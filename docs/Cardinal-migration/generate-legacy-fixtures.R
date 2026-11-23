# Rebuild compact scientific references using the preserved pre-migration function.
# Run from the repository root only when deliberately rebuilding the reference:
#   Rscript --vanilla docs/Cardinal-migration/generate-legacy-fixtures.R
#
# This script OVERWRITES the saved fixture. Ordinary validation reads that fixture
# instead; never rerun this script to hide a numerical regression. The frozen
# source is unchanged, and incidental random/time-based metadata IDs are omitted.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
legacy <- new.env(parent = asNamespace("moleculaR"))
# Bind only the names needed by the unchanged historical source. These bindings
# stay inside this reference-generation environment, never the package namespace.
legacy$ppp <- spatstat.geom::ppp
legacy$owin <- spatstat.geom::owin
legacy$as.polygonal <- spatstat.geom::as.polygonal
sys.source("tests/fixtures/legacy/searchAnalyte.R", envir = legacy)

# An asymmetric tissue domain with a hole catches flips, coordinate reordering,
# and accidental loss of acquired pixels whose analyte intensity is zero.
coords <- expand.grid(x = 1:16, y = 1:11)
coords <- coords[!(coords$x <= 2 & coords$y >= 9) & !(coords$x == 12 & coords$y == 3), ,
    drop = FALSE]
rownames(coords) <- NULL
n <- nrow(coords)
mz <- c(99, 99.9, 100, 100.1, 101)

# Place a bright central region in the analyte signal. Omit selected neighboring
# peaks in some pixels to exercise both Gaussian weighting and mean denominators.
signal <- ifelse(coords$x >= 6 & coords$x <= 11 & coords$y >= 4 & coords$y <= 8, 20, 1) +
    (seq_len(n)%%5)/10
intensity <- rbind(rep(2, n), signal * 0.1, signal, signal * 0.2, rep(3, n))
intensity[2, seq(1, n, 3)] <- 0
intensity[4, seq(1, n, 4)] <- 0
intensity[2:4, n] <- 0

# Build the historical sparse representation only for reference calculation.
spectra <- lapply(seq_len(n), function(i) {
    keep <- intensity[, i] > 0
    MALDIquant::createMassPeaks(mass = mz[keep], intensity = intensity[keep, i], metaData = list(imaging = list(pos = unlist(coords[i,
        ]))))
})
spData <- createSparseMat(spectra)
spwin <- createSpatialWindow(coords, plot = FALSE)
patterns <- lapply(c("Gaussian", "sum", "mean", "max"), function(method) legacy$searchAnalyte(100,
    0.1, spData, wMethod = method, spwin = spwin, verifiedMasses = c(100, 200), metaData = list(label = "fixture")))
names(patterns) <- c("Gaussian", "sum", "mean", "max")

# Fix both bandwidth and random seed so density and significance comparisons
# isolate spectral handling from spatial stochastic or bandwidth-selection effects.
mpm <- suppressWarnings(probMap(patterns$Gaussian, bw = 1, seed = 1, verbose = FALSE))
summaryMpm <- function(x) list(bw = x$bw, rhoMoi = x$rhoMoi$v, rhoCsr = x$rhoCsr$v, hotspot = x$hotspotIm$v,
    coldspot = x$coldspotIm$v, hotspotCoordinates = cbind(x$hotspotpp$x, x$hotspotpp$y),
    coldspotCoordinates = cbind(x$coldspotpp$x, x$coldspotpp$y))
reference <- list(startingCommit = "ce10618b608e840f888bb71bb68794c8b44535ab", sourceMd5 = unname(tools::md5sum("tests/fixtures/legacy/searchAnalyte.R")),
    mz = mz, intensity = intensity, coordinates = coords, target = 100, fwhm = 0.1, patterns = lapply(patterns,
        function(x) list(x = x$x, y = x$y, intensity = x$marks$intensity, mzConfirmed = x$metaData$mzConfirmed)),
    mpm = summaryMpm(mpm), defaultFwhm = getFwhm(useDefaultFwhm(), c(249, 400, 788.5447,
        1000)), versions = vapply(c("MALDIquant", "spatstat.geom", "spatstat.explore"), function(x) as.character(packageVersion(x)),
        character(1)))

saveRDS(reference, "tests/fixtures/legacy/single-analyte.rds", version = 3)
cat("Legacy reference saved:", n, "pixels;", patterns$Gaussian$n, "detections\n")
print(reference$versions)
cat("ProtGenerics processingData exported:", "processingData" %in% getNamespaceExports("ProtGenerics"),
    "\n")
