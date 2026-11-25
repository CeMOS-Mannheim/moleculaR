# Deliberately capture the Phase 3 baseline before modifying preprocessing.
# Run from the repository root with the isolated migration library available.
# This OVERWRITES phase4 fixtures and the extracted numeric FTICR reference;
# ordinary regression tests must only read these artifacts.
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
baseline <- "2308bf4"
destination <- "tests/fixtures/legacy/phase4"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
sourceFiles <- c("estimateFwhm.R", "normalizeIntensity.R", "filterPeaks.R",
                 "readSingleSpect.R", "useDefaultFwhm.R")
legacy <- new.env(parent = asNamespace("moleculaR"))
for (file in sourceFiles) {
    source <- system2("git", c("show", paste0(baseline, ":R/", file)), stdout = TRUE)
    stopifnot(is.null(attr(source, "status")))
    target <- file.path(destination, file)
    writeLines(source, target)
    sys.source(target, legacy)
}

example <- new.env(parent = baseenv())
utils::data("processed-example-Data", package = "moleculaR", envir = example)
model <- example$fwhmObj
interpolation <- mget(c("x", "y", "method", "yleft", "yright", "f", "na.rm"),
                      envir = environment(model$fwhmInterpolator))
reference <- list(sourceCommit = baseline,
    sourceMd5 = unname(tools::md5sum("data/processed-example-Data.RData")),
    peaks = model$peaks, measuredFwhm = model$fwhmVals, interpolation = interpolation)
saveRDS(reference, "inst/extdata/fticr-fwhm-reference.rds", version = 3)
referenceMasses <- c(50, 100, 150, 249, 400, 788.5447, 1000, 1400, 2000)

axis <- c(99.9, 100, 100.1, 200, 300)
signal <- cbind(c(2, 10, 4, 5, 0), c(4, 20, 8, 0, 6), c(1, 5, 2, 10, 12))
coordinates <- data.frame(x = 1:3, y = 1)
makeLegacy <- function(dropZeros) lapply(seq_len(ncol(signal)), function(i) {
    keep <- if (dropZeros) signal[, i] > 0 else rep(TRUE, length(axis))
    MALDIquant::createMassPeaks(mass = axis[keep], intensity = signal[keep, i],
        metaData = list(imaging = list(pos = unlist(coordinates[i, ]))))
})
spectra <- makeLegacy(TRUE)
normalizations <- lapply(c("TIC", "RMS", "median", "IS", "medianFoldChange"), function(method) {
    normalized <- legacy$normalizeIntensity(spectra, method, mz = 200, tol = 0)
    lapply(normalized, MALDIquant::intensity)
})
names(normalizations) <- c("TIC", "RMS", "median", "IS", "medianFoldChange")
fold <- lapply(c("mean", "median", "sum"), function(method) {
    lapply(legacy$normalizeIntensity(spectra, "medianFoldChange", meanMethod = method),
           MALDIquant::intensity)
})
names(fold) <- c("mean", "median", "sum")
fullFold <- lapply(legacy$normalizeIntensity(makeLegacy(FALSE), "medianFoldChange"),
                   MALDIquant::intensity)
filtered <- lapply(legacy$filterPeaks(spectra, minFreq = 2/3), MALDIquant::mass)

mass <- seq(100, 200, by = 0.002)
centres <- seq(110, 190, by = 4)
intensity <- 1 + Reduce(`+`, lapply(centres, function(m) 100 * exp(-0.5 * ((mass - m)/0.05)^2)))
profile <- MALDIquant::createMassSpectrum(mass, intensity)
peaks <- MALDIquant::detectPeaks(profile, method = "SuperSmoother", SNR = 5)
widths <- legacy$.getFwhm(profile, peaks)
set.seed(4)
fit <- legacy$estimateFwhm(profile)

# Capture representative real-profile detections and fitted widths without
# copying the large input TSV into the fixture.
profilePath <- "inst/extdata/pos-XIII-82492-singleSpectrum.tsv"
real <- legacy$readSingleSpect(profilePath)
realPeaks <- MALDIquant::detectPeaks(real, method = "SuperSmoother", SNR = 5)
realWidths <- legacy$.getFwhm(real, realPeaks)
set.seed(4)
realFit <- legacy$estimateFwhm(real)
realMasses <- c(150, 249, 400, 788.5447, 1000)

fixture <- list(sourceCommit = baseline,
    sourceMd5 = tools::md5sum(file.path(destination, sourceFiles)),
    axis = axis, signal = signal, coordinates = coordinates,
    normalizations = normalizations, fold = fold, fullFold = fullFold, filtered = filtered,
    fticrMasses = referenceMasses,
    fticrPredictions = getFwhm(legacy$useDefaultFwhm(), referenceMasses),
    profile = list(mz = mass, intensity = intensity),
    peakIndices = MALDIquant::match.closest(MALDIquant::mass(peaks), mass),
    widths = widths, predictionMasses = centres, predictions = getFwhm(fit, centres),
    real = list(inputMd5 = unname(tools::md5sum(profilePath)),
        peakMz = MALDIquant::mass(realPeaks), widths = realWidths,
        predictionMasses = realMasses, predictions = getFwhm(realFit, realMasses)))
saveRDS(fixture, file.path(destination, "preprocessing.rds"), version = 3)
cat("Saved immutable phase4 baseline:", nrow(widths), "synthetic peaks,",
    nrow(realWidths), "real peaks; native FTICR reference", file.size("inst/extdata/fticr-fwhm-reference.rds"), "bytes\n")
