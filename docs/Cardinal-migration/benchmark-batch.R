# Compare Phase 3 with repeated single extraction in fresh R processes.
# Run after installing the current checkout into the local migration library:
#   Rscript --vanilla docs/Cardinal-migration/benchmark-batch.R [imzML]
# Input files are read only. Reports describe one local dataset, uncontrolled
# disk-cache conditions and whole-process peak working set, not incremental RAM.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
stopifnot(requireNamespace("Matrix", quietly = TRUE))
root <- normalizePath(getwd(), winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
path <- normalizePath(if (length(args)) args[1] else
    "data-external/pos-XIII-82492-msi.imzML", winslash = "/", mustWork = TRUE)
masses <- c(734.5694, 760.5851, 782.5673, 784.5851, 788.5447, 810.5290, 832.5824, 834.6007)

measure <- function(root, path, mode, pixels, chunkSize, workers, masses) {
    suppressPackageStartupMessages(library(moleculaR))
    imported <- system.time(x <- suppressWarnings(readMsi(path, memory = FALSE, verbose = FALSE)))
    n <- moleculaR:::.msiNumberOfPixels(x)
    if (!is.na(pixels)) x <- moleculaR:::.msiSelectPixels(x, seq_len(min(n, pixels)))
    n <- moleculaR:::.msiNumberOfPixels(x)
    widths <- getFwhm(useDefaultFwhm(), masses)
    backend <- moleculaR:::.msiBatchBackend(workers, NULL)
    cache <- if (mode == "cache-hit") new.env(parent = emptyenv()) else NULL
    if (!is.null(cache)) {
        invisible(moleculaR:::.msiExtractIons(x, masses, widths, cache = cache))
    }
    allocationLog <- tempfile()
    Rprofmem(allocationLog, threshold = 1e6)
    timing <- system.time({
        values <- if (mode == "single") {
            Matrix::Matrix(vapply(seq_along(masses), function(i) {
                moleculaR:::.msiExtractIon(x, masses[i], widths[i])
            }, numeric(n)), sparse = TRUE)
        } else {
            moleculaR:::.msiExtractIons(x, masses, widths, chunkSize = chunkSize,
                BPPARAM = backend, cache = cache)
        }
    })
    Rprofmem(NULL)
    allocations <- suppressWarnings(as.numeric(sub(" .*", "", readLines(allocationLog))))
    unlink(allocationLog)
    peak <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
        c("-NoProfile", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64",
                                                  Sys.getpid()))), stdout = TRUE)) else NA_real_
    report <- data.frame(mode = mode, pixels = n, targets = length(masses),
        chunkSize = if (is.null(chunkSize)) NA_real_ else chunkSize, workers = workers,
        importSeconds = unname(imported["elapsed"]), extractionSeconds = unname(timing["elapsed"]),
        resultBytes = as.numeric(object.size(values)), detections = length(values@x),
        peakWorkingSetBytes = peak,
        largestRAllocationBytes = if (any(is.finite(allocations))) max(allocations, na.rm = TRUE) else 0,
        ibdBytes = file.size(paste0(tools::file_path_sans_ext(path), ".ibd")),
        R = getRversion() |> as.character(), Cardinal = as.character(packageVersion("Cardinal")),
        matter = as.character(packageVersion("matter")), platform = R.version$platform)
    list(report = report, values = values)
}

variants <- list(
    list(mode = "single", pixels = NA_real_, chunkSize = NULL, workers = 1),
    list(mode = "batch", pixels = NA_real_, chunkSize = NULL, workers = 1),
    list(mode = "cache-hit", pixels = NA_real_, chunkSize = NULL, workers = 1),
    list(mode = "batch", pixels = 2000, chunkSize = 100, workers = 1),
    list(mode = "batch", pixels = 2000, chunkSize = 500, workers = 1),
    list(mode = "batch", pixels = 2000, chunkSize = 2000, workers = 1),
    list(mode = "batch", pixels = 2000, chunkSize = 2000, workers = 2)
)
reports <- vector("list", length(variants))
artifacts <- file.path(root, "renv/local/phase3-benchmark")
dir.create(artifacts, recursive = TRUE, showWarnings = FALSE)
reference <- NULL
for (i in seq_along(variants)) {
    variant <- variants[[i]]
    cat("Measuring", variant$mode, "pixels", variant$pixels, "workers", variant$workers, "\n")
    flush.console()
    result <- callr::r(measure, c(list(root = root, path = path), variant, list(masses = masses)),
        libpath = .libPaths(), env = c(LC_ALL = "C", RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE"))
    saveRDS(result, file.path(artifacts, paste0("variant-", i, ".rds")))
    if (is.null(reference)) reference <- result$values
    expected <- reference[seq_len(nrow(result$values)), , drop = FALSE]
    stopifnot(isTRUE(all.equal(result$values, expected, tolerance = 1e-10)))
    reports[[i]] <- result$report
    print(result$report[, c("mode", "pixels", "chunkSize", "workers", "extractionSeconds", "peakWorkingSetBytes")])
}
write.csv(do.call(rbind, reports), "docs/Cardinal-migration/batch-benchmark.csv", row.names = FALSE)
