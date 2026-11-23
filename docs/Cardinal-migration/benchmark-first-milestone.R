# Benchmark an existing file-backed MSI dataset in a fresh R process.
# Run from the repository root, with an optional imzML path:
#   Rscript --vanilla docs/Cardinal-migration/benchmark-first-milestone.R [path]
#
# The input pair is read only. CSV timing/memory observations go beside this
# script; allocation traces stay in ignored renv/local storage. These timings
# do not establish large-scale or cross-platform performance guarantees.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "data-external/pos-XIII-82492-msi.imzML"
stopifnot(file.exists(path))

# Measure import and a first Gaussian extraction separately. No global sparse
# representation is created, and physical coordinate scaling is not applied.
import <- system.time(x <- readMsi(path, memory = FALSE, verbose = FALSE))
coordinates <- moleculaR:::.msiCoordinates(x)
extract <- system.time(values <- moleculaR:::.msiExtractIon(x, 788.5447, getFwhm(useDefaultFwhm(),
    788.5447)))
invisible(gc(reset = TRUE))

# Repeat with allocation profiling. The threshold reports individual R
# allocations >= 1 MB, not total allocations or native-library memory use.
allocationLog <- "renv/local/cardinal-extraction-allocations.log"
dir.create(dirname(allocationLog), recursive = TRUE, showWarnings = FALSE)
Rprofmem(allocationLog, threshold = 1e+06)
repeatTiming <- system.time(repeated <- moleculaR:::.msiExtractIon(x, 788.5447, getFwhm(useDefaultFwhm(),
    788.5447)))
Rprofmem(NULL)
stopifnot(identical(values, repeated))

# Read the peak working set of this actual R process, not the Rscript launcher.
# It includes import/package-loading overhead as well as both extractions.
peakBytes <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe", c("-NoProfile",
    "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))),
    stdout = TRUE)) else NA_real_
allocations <- readLines(allocationLog)
allocationBytes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
maxAllocation <- if (any(is.finite(allocationBytes))) max(allocationBytes, na.rm = TRUE) else 0

# Retain raw sizes and units in the CSV so reports can distinguish process
# memory, pointer-object size, disk storage, and extraction time.
report <- data.frame(dataset = basename(path), pixels = nrow(coordinates), imzmlBytes = file.size(path),
    ibdBytes = file.size(paste0(tools::file_path_sans_ext(path), ".ibd")), importSeconds = unname(import["elapsed"]),
    extractionSeconds = unname(extract["elapsed"]), repeatSeconds = unname(repeatTiming["elapsed"]),
    detectedPixels = sum(values > 0), backed = methods::is(Cardinal::intensity(x), "matter"),
    objectBytes = as.numeric(object.size(x)), peakWorkingSetBytes = peakBytes, largestProfiledAllocationBytes = maxAllocation)
print(report)
write.csv(report, "docs/Cardinal-migration/first-milestone-benchmark.csv", row.names = FALSE)
