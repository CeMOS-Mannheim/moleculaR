# Inspect Cardinal methods and exercise a tiny native MSI object.
# Run from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/inspect-api.R
#
# This diagnostic prints method bodies to help verify version-specific behavior.
# It does not alter package code, input files, or the dependency lockfile.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
suppressPackageStartupMessages(library(Cardinal))

# Inspect dispatch for both supported representations before relying on an API.
print(grep("process|Process|centroid|run|coord", getNamespaceExports("Cardinal"), value = TRUE))
for (fun in c("spectrapply", "process", "centroided", "mz", "intensity", "coord", "run")) {
    cat("\n", fun, "\n")
    for (cls in c("MSImagingArrays", "MSImagingExperiment")) print(try(methods::selectMethod(fun,
        cls), silent = TRUE))
}
print(Cardinal::readImzML)
print(Cardinal::PositionDataFrame)

# Two spectra with different mass axes demonstrate acquisition-order access.
x <- Cardinal::MSImagingArrays(spectraData = list(mz = list(c(100, 200), c(100, 201)), intensity = list(c(1,
    2), c(3, 4))), pixelData = Cardinal::PositionDataFrame(data.frame(x = 1:2, y = c(1, 1))),
    centroided = TRUE, continuous = FALSE)
print(x)
print(methods::slotNames(x))
print(Cardinal::coord(x))
print(Cardinal::run(x))

# Explicit serial iteration returns one aggregate per spectrum. Normalization
# below is queued only; the printed processing state illustrates that distinction.
print(Cardinal::spectrapply(x, function(i, mz) sum(i), index = "mz", BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE))
print(Cardinal::normalize(x, method = "tic"))

# Writer precision matters when constructing numerical round-trip fixtures.
cat("\nCardinalIO writer numeric types:\n")
print(methods::selectMethod("writeImzML", "ImzML"))
