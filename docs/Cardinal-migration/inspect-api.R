.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
suppressPackageStartupMessages(library(Cardinal))
print(grep("process|Process|centroid|run|coord", getNamespaceExports("Cardinal"), value = TRUE))
for (fun in c("spectrapply", "process", "centroided", "mz", "intensity", "coord", "run")) {
    cat("\n", fun, "\n")
    for (cls in c("MSImagingArrays", "MSImagingExperiment"))
        print(try(methods::selectMethod(fun, cls), silent = TRUE))
}
print(Cardinal::readImzML)
print(Cardinal::PositionDataFrame)
x <- Cardinal::MSImagingArrays(
    spectraData = list(mz = list(c(100, 200), c(100, 201)),
                       intensity = list(c(1, 2), c(3, 4))),
    pixelData = Cardinal::PositionDataFrame(data.frame(x = 1:2, y = c(1, 1))),
    centroided = TRUE, continuous = FALSE)
print(x)
print(methods::slotNames(x))
print(Cardinal::coord(x))
print(Cardinal::run(x))
print(Cardinal::spectrapply(x, function(i, mz) sum(i), index = "mz",
                          BPPARAM = BiocParallel::SerialParam(), verbose = FALSE))
print(Cardinal::normalize(x, method = "tic"))
cat("\nCardinalIO writer numeric types:\n")
print(methods::selectMethod("writeImzML", "ImzML"))
