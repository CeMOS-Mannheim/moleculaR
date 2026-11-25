# Quantify native preprocessing differences against the frozen Phase 3 baseline.
# Run from the repository root; raw spectra are read only. CSVs record changes
# rather than asserting equivalence between different automatic peak pickers.
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
pkgload::load_all(".", quiet = TRUE)
f <- readRDS("tests/fixtures/legacy/phase4/preprocessing.rds")
profilePath <- "inst/extdata/pos-XIII-82492-singleSpectrum.tsv"
profile <- readSingleSpect(profilePath)
backend <- BiocParallel::SerialParam()
indices <- moleculaR:::.msiPickProfile(profile, "diff", 5, backend)
newMz <- profile$mz[indices]
newWidths <- moleculaR:::.estimatePeakWidths(profile$mz, profile$intensity, indices)
common <- match(newMz, f$real$peakMz)
matched <- !is.na(common)
widthDifference <- newWidths[matched] - f$real$widths$fwhmValues[common[matched]]
set.seed(4)
timing <- system.time(model <- estimateFwhm(profile))
newPredictions <- getFwhm(model, f$real$predictionMasses)
comparison <- data.frame(mz = f$real$predictionMasses, legacyWidth = f$real$predictions,
    CardinalWidth = newPredictions, absoluteDifference = newPredictions - f$real$predictions,
    relativeDifference = (newPredictions - f$real$predictions)/f$real$predictions)
write.csv(comparison, "docs/Cardinal-migration/phase4-fwhm-comparison.csv", row.names = FALSE)
summary <- data.frame(profileSource = basename(profilePath),
    profileMd5 = unname(tools::md5sum(profilePath)),
    legacyPeaks = length(f$real$peakMz), CardinalPeaks = length(newMz),
    commonPeaks = sum(matched), legacyOnly = sum(!f$real$peakMz %in% newMz),
    CardinalOnly = sum(!matched),
    maximumCommonWidthDifference = max(abs(widthDifference), na.rm = TRUE),
    retainedFitWidths = length(model$peaks), discardedWidths = sum(model$estimation$discardedWidths),
    fitSeconds = unname(timing["elapsed"]))
write.csv(summary, "docs/Cardinal-migration/phase4-peak-comparison.csv", row.names = FALSE)
print(summary)
print(comparison)

path <- "data-external/pos-XIII-82492-msi.imzML"
if (file.exists(path)) {
    x <- suppressWarnings(readMsi(path, verbose = FALSE))
    target <- 788.5447
    oldWidth <- f$real$predictions[match(target, f$real$predictionMasses)]
    newWidth <- getFwhm(model, target)
    old <- moleculaR:::.msiExtractIon(x, target, oldWidth)
    new <- moleculaR:::.msiExtractIon(x, target, newWidth)
    domain <- moleculaR:::.msiBatchDomain(x, NA)
    pointPattern <- function(values) {
        keep <- which(values > 0)
        analytePointPattern(x = domain$coordinates$x[keep], y = domain$coordinates$y[keep],
            win = domain$window, intensity = values[keep], mzVals = target)
    }
    first <- suppressWarnings(probMap(pointPattern(old), bw = 2, seed = 1, verbose = FALSE))
    second <- suppressWarnings(probMap(pointPattern(new), bw = 2, seed = 1, verbose = FALSE))
    classify <- function(mpm) ifelse(mpm$hotspotIm$v > 0, 1L, ifelse(mpm$coldspotIm$v > 0, -1L, 0L))
    maps <- data.frame(msiSource = basename(path), msiRepresentation = "provided reduced centroided MSI",
        widthSource = basename(profilePath), mz = target, legacyWidth = oldWidth, CardinalWidth = newWidth,
        acquiredPixels = length(old), legacyDetections = sum(old > 0), CardinalDetections = sum(new > 0),
        changedIntensities = sum(old != new), maxIntensityDifference = max(abs(old - new)),
        maxDensityDifference = max(abs(first$rhoMoi$v - second$rhoMoi$v), na.rm = TRUE),
        changedClassifications = sum(classify(first) != classify(second), na.rm = TRUE),
        bandwidth = 2, seed = 1)
    write.csv(maps, "docs/Cardinal-migration/phase4-map-comparison.csv", row.names = FALSE)
    print(maps)

    # A moderate subset of the provided reduced centroided data demonstrates
    # streamed normalization. It is never used for picking or width estimation.
    subset <- moleculaR:::.msiSelectPixels(x, seq_len(2000))
    allocations <- tempfile()
    Rprofmem(allocations, threshold = 1e6)
    elapsed <- system.time(normalized <- normalizeIntensity(subset, "TIC"))
    Rprofmem(NULL)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(allocations))))
    sources <- moleculaR:::.msiBackingState(normalized)
    peak <- if (.Platform$OS.type == "windows") as.numeric(system2("powershell.exe",
        c("-NoProfile", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid()))),
        stdout = TRUE)) else NA_real_
    benchmark <- data.frame(msiSource = basename(path), msiRepresentation = "provided reduced centroided MSI",
        normalizationBasis = "available reduced feature intensities", pixels = 2000,
        method = "Cardinal TIC", elapsedSeconds = unname(elapsed["elapsed"]),
        fileBacked = methods::is(Cardinal::intensity(normalized), "matter"),
        outputBinaryBytes = sum(sources$size[!sources$path %in% moleculaR:::.msiBackingState(subset)$path]),
        objectBytes = as.numeric(object.size(normalized)),
        largestRAllocationBytes = if (any(is.finite(sizes))) max(sizes, na.rm = TRUE) else 0,
        wholeProcessPeakBytes = peak, Cardinal = as.character(packageVersion("Cardinal")),
        matter = as.character(packageVersion("matter")))
    write.csv(benchmark, "docs/Cardinal-migration/phase4-benchmark.csv", row.names = FALSE)
    print(benchmark)
}
