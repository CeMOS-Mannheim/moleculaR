# Synthetic MSI fixtures deliberately retain an asymmetric pixel order. This
# catches accidental coordinate sorting as well as mass-by-pixel transposition.
# None of these constructors uses MALDIquant or a global sparse conversion.
make_msi_arrays <- function(mz = list(c(99.9, 100, 100.1), c(100, 101), numeric()), intensity = list(c(2,
    10, 4), c(3, 9), numeric()), coordinates = data.frame(x = c(3, 1, 2), y = c(1, 2, 1)),
    runs = rep("sample", nrow(coordinates)), centroided = TRUE) {
    Cardinal::MSImagingArrays(spectraData = list(mz = mz, intensity = intensity), pixelData = Cardinal::PositionDataFrame(coordinates,
        run = runs), centroided = centroided, continuous = FALSE)
}

make_msi_experiment <- function(mz = c(99.9, 100, 100.1, 101), intensity = cbind(c(2, 10,
    4, 0), c(0, 3, 0, 9), c(0, 0, 0, 0)), coordinates = data.frame(x = c(3, 1, 2), y = c(1,
    2, 1))) {
    Cardinal::MSImagingExperiment(spectraData = list(intensity = intensity), featureData = Cardinal::MassDataFrame(mz),
        pixelData = Cardinal::PositionDataFrame(coordinates), centroided = TRUE)
}

# The saved values are an immutable scientific reference from the starting
# implementation. Do not regenerate them simply to make a regression pass.
legacy_fixture <- function() readRDS(testthat::test_path("..", "fixtures", "legacy", "single-analyte.rds"))

# Replay the frozen source without changing it. Its historical unqualified
# geometry calls are bound only in this private test environment; production
# functions use explicit namespaces and have no dynamic importer.
legacy_search_environment <- function() {
    legacy <- new.env(parent = asNamespace("moleculaR"))
    legacy$ppp <- spatstat.geom::ppp
    legacy$owin <- spatstat.geom::owin
    legacy$as.polygonal <- spatstat.geom::as.polygonal

    sys.source(testthat::test_path("..", "fixtures", "legacy", "searchAnalyte.R"), envir = legacy)

    legacy
}

fixture_cardinal <- function(arrays = TRUE, dropZeros = FALSE) {
    f <- legacy_fixture()
    if (!arrays)
        return(make_msi_experiment(f$mz, f$intensity, f$coordinates))
    mz <- lapply(seq_len(ncol(f$intensity)), function(i) f$mz[if (dropZeros)
        f$intensity[, i] > 0 else rep(TRUE, length(f$mz))])
    intensity <- lapply(seq_len(ncol(f$intensity)), function(i) f$intensity[if (dropZeros)
        f$intensity[, i] > 0 else rep(TRUE, length(f$mz)), i])
    make_msi_arrays(mz, intensity, f$coordinates)
}

# Preserve the double-precision fixture through disk serialization. CardinalIO
# otherwise defaults to float32 intensities and changes the numerical input.
write_msi_pair <- function(x, path) {
    Cardinal::writeMSIData(x, file = path, bundle = FALSE, verbose = FALSE, intensity.type = "float64",
        BPPARAM = BiocParallel::SerialParam())
    path
}
