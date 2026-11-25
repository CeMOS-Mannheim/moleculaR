# Exercise the same session data/analysis operations used by the Shiny modules.
# Rscript --vanilla docs/Cardinal-migration/validate-phase5.R [optional imzML]
# The optional acquisition must be known-centroided. It is never peak-picked
# or used for measured widths; the approximate FTICR model is explicit.
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
example <- normalizePath("inst/extdata/shiny-example", winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
raw <- if (length(args)) args[1] else "data-external/pos-XIII-82492-msi.imzML"

evaluate <- function(source, path = NULL) {
    state <- new.env(parent = emptyenv())
    state$root <- moleculaR:::.shinySessionDirectory()
    state$auxiliary <- list()
    on.exit(moleculaR:::.shinyRemove(state$root, state$root, TRUE))
    paths <- if (source == "upload") c(path, sub("[.]imzML$", ".ibd", path)) else file.path(example, c("example.imzML", "example.ibd"))
    before <- tools::md5sum(paths)
    input <- list(imzmlFile = if (source == "upload") data.frame(name = basename(paths), datapath = paths) else NULL,
        fwhmSource = "default", fwhmModality = "FTICR", normalization = "none",
        mz = 788.5447, mzStr = "760.5851,788.5447", wMethod = "Gaussian", confirmedOnly = source == "example",
        standardize = FALSE, seed = 1, bandwidthMode = "fixed", bandwidth = 2,
        lipidView = "class", lipidClass = "PC(x:x)", adducts = "[M+H]+")
    initialization <- system.time(moleculaR:::.shinyInitialize(state, input, source, example))
    log <- tempfile(tmpdir = "renv/local", fileext = ".log")
    Rprofmem(log, threshold = 1e6)
    singleTime <- system.time(single <- suppressWarnings(moleculaR:::.shinyAnalyze(state, input, "single")))
    Rprofmem(NULL)
    batchTime <- system.time(batch <- suppressWarnings(moleculaR:::.shinyAnalyze(state, input, "batch")))
    lipidTime <- if (source == "example") system.time(lipid <- suppressWarnings(moleculaR:::.shinyAnalyze(state, input, "lipid"))) else c(elapsed = NA_real_)
    reference <- suppressWarnings(probMap(batchAnalyteSearch(state$msData, state$fwhmObj, state$spwin,
        input$mz, verifiedMasses = state$verifiedMasses, confirmedOnly = input$confirmedOnly, verbose = FALSE),
        bw = 2, seed = 1, verbose = FALSE))
    differences <- vapply(c("rhoMoi", "rhoCsr", "hotspotIm", "coldspotIm"), function(field) {
        max(abs(single$map[[field]]$v - reference[[field]]$v), na.rm = TRUE)
    }, numeric(1))
    stopifnot(all(differences == 0), methods::is(Cardinal::intensity(state$msData), "matter"),
        identical(before, tools::md5sum(paths)))
    prefix <- paste0("renv/local/phase5-", source)
    for (format in c("rds", "csv", "png")) moleculaR:::.shinySaveAnalysis(single, paste0(prefix, ".", format), format)
    allocations <- suppressWarnings(as.numeric(sub(" .*", "", readLines(log))))
    peak <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe", c("-NoProfile", "-Command",
        shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))), stdout = TRUE)) else NA_real_
    data.frame(source = source, filename = basename(paths[1]), centroidedReduced = TRUE,
        pixels = moleculaR:::.msiNumberOfPixels(state$msData), imzmlBytes = file.size(paths[1]), ibdBytes = file.size(paths[2]),
        initializeSeconds = unname(initialization["elapsed"]), singleSeconds = unname(singleTime["elapsed"]),
        batchSeconds = unname(batchTime["elapsed"]), lipidSeconds = unname(lipidTime["elapsed"]),
        singlePoints = single$pointPattern$n, batchAnalytes = nrow(batch$pointPattern$metaData),
        lipidAnalytes = if (source == "example") nrow(lipid$pointPattern$metaData) else NA_integer_,
        backed = methods::is(Cardinal::intensity(state$msData), "matter"), objectBytes = as.numeric(object.size(state$msData)),
        largestProfiledRAllocationBytes = if (any(is.finite(allocations))) max(allocations, na.rm = TRUE) else 0,
        peakProcessWorkingSetBytes = peak, maxMapDifference = max(differences), originalFilesUnchanged = TRUE)
}

report <- evaluate("example")
if (file.exists(raw)) report <- rbind(report, evaluate("upload", normalizePath(raw, winslash = "/")))
print(report)
utils::write.csv(report, "docs/Cardinal-migration/phase5-benchmark.csv", row.names = FALSE)
