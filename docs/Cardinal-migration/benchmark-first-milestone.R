# Fresh-process benchmark of an existing local file-backed MSI dataset.
# Invoke with optional path; raw user datasets are never changed or checked in.
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "data-external/pos-XIII-82492-msi.imzML"
stopifnot(file.exists(path))
import <- system.time(x <- readMsi(path, memory = FALSE, verbose = FALSE))
coordinates <- moleculaR:::.msiCoordinates(x)
extract <- system.time(values <- moleculaR:::.msiExtractIon(x, 788.5447,
    getFwhm(useDefaultFwhm(), 788.5447)))
invisible(gc(reset = TRUE))
allocationLog <- "renv/local/cardinal-extraction-allocations.log"
dir.create(dirname(allocationLog), recursive = TRUE, showWarnings = FALSE)
Rprofmem(allocationLog, threshold = 1e6)
repeatTiming <- system.time(repeated <- moleculaR:::.msiExtractIon(x, 788.5447,
    getFwhm(useDefaultFwhm(), 788.5447)))
Rprofmem(NULL)
stopifnot(identical(values, repeated))
peakBytes <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
    c("-NoProfile", "-Command", shQuote(sprintf(
        "(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))), stdout = TRUE)) else NA_real_
allocations <- readLines(allocationLog)
allocationBytes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
maxAllocation <- if (any(is.finite(allocationBytes))) max(allocationBytes, na.rm = TRUE) else 0
report <- data.frame(dataset = basename(path), pixels = nrow(coordinates),
    imzmlBytes = file.size(path),
    ibdBytes = file.size(paste0(tools::file_path_sans_ext(path), ".ibd")),
    importSeconds = unname(import["elapsed"]),
    extractionSeconds = unname(extract["elapsed"]),
    repeatSeconds = unname(repeatTiming["elapsed"]),
    detectedPixels = sum(values > 0),
    backed = methods::is(Cardinal::intensity(x), "matter"),
    objectBytes = as.numeric(object.size(x)), peakWorkingSetBytes = peakBytes,
    largestProfiledAllocationBytes = maxAllocation)
print(report)
write.csv(report, "docs/Cardinal-migration/first-milestone-benchmark.csv", row.names = FALSE)
