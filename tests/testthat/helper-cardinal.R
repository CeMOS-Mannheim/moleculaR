make_msi_arrays <- function(mz = list(c(99.9, 100, 100.1), c(100, 101), numeric()),
                            intensity = list(c(2, 10, 4), c(3, 9), numeric()),
                            coordinates = data.frame(x = c(3, 1, 2), y = c(1, 2, 1)),
                            runs = rep("sample", nrow(coordinates)), centroided = TRUE) {
    Cardinal::MSImagingArrays(spectraData = list(mz = mz, intensity = intensity),
        pixelData = Cardinal::PositionDataFrame(coordinates, run = runs),
        centroided = centroided, continuous = FALSE)
}

make_msi_experiment <- function(mz = c(99.9, 100, 100.1, 101),
                                intensity = matrix(c(2, 10, 4, 0, 0, 3, 0, 9, 0, 0, 0, 0), 4),
                                coordinates = data.frame(x = c(3, 1, 2), y = c(1, 2, 1))) {
    Cardinal::MSImagingExperiment(spectraData = list(intensity = intensity),
        featureData = Cardinal::MassDataFrame(mz),
        pixelData = Cardinal::PositionDataFrame(coordinates), centroided = TRUE)
}

legacy_fixture <- function() readRDS(testthat::test_path("..", "fixtures", "legacy", "single-analyte.rds"))

fixture_cardinal <- function(arrays = TRUE, dropZeros = FALSE) {
    f <- legacy_fixture()
    if (!arrays) return(make_msi_experiment(f$mz, f$intensity, f$coordinates))
    mz <- lapply(seq_len(ncol(f$intensity)), function(i)
        f$mz[if (dropZeros) f$intensity[, i] > 0 else rep(TRUE, length(f$mz))])
    intensity <- lapply(seq_len(ncol(f$intensity)), function(i)
        f$intensity[if (dropZeros) f$intensity[, i] > 0 else rep(TRUE, length(f$mz)), i])
    make_msi_arrays(mz, intensity, f$coordinates)
}

write_msi_pair <- function(x, path) {
    Cardinal::writeMSIData(x, file = path, bundle = FALSE, verbose = FALSE,
                          intensity.type = "float64", BPPARAM = BiocParallel::SerialParam())
    path
}
