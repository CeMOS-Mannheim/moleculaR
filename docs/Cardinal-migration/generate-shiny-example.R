# Deliberately prepare a small native example from the preserved reduced data.
# Development only: runtime apps read this imzML pair and never load legacy
# spectra. Run from the repository root; this overwrites the example artifacts.
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
legacy <- new.env(parent = baseenv())
utils::data("processed-example-Data", package = "moleculaR", envir = legacy)
coordinates <- as.data.frame(MALDIquant::coordinates(legacy$msData))
centre <- vapply(coordinates, stats::median, numeric(1))
keep <- which(coordinates$x >= floor(centre[1]) - 15 & coordinates$x <= floor(centre[1]) + 16 &
              coordinates$y >= floor(centre[2]) - 15 & coordinates$y <= floor(centre[2]) + 16)
stopifnot(length(keep) >= 400)
native <- Cardinal::MSImagingArrays(
    spectraData = list(mz = lapply(legacy$msData[keep], MALDIquant::mass),
                       intensity = lapply(legacy$msData[keep], MALDIquant::intensity)),
    pixelData = Cardinal::PositionDataFrame(coordinates[keep, ], run = "glioblastoma-example"),
    centroided = TRUE, continuous = FALSE)
output <- "inst/extdata/shiny-example"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
Cardinal::writeMSIData(native, file = file.path(output, "example.imzML"), bundle = FALSE,
    intensity.type = "float64", verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
annotations <- data.frame(mz = as.numeric(legacy$mtspc$mz), fdr = as.numeric(legacy$mtspc$fdr))
utils::write.csv(annotations, file.path(output, "annotations.csv"), row.names = FALSE)
model <- useDefaultFwhm()
database <- legacy$sldb
eligible <- rep(FALSE, nrow(database))
for (column in c("Exact m/z of [M+H]+", "Exact m/z of [M+Na]+", "Exact m/z of [M+K]+")) {
    masses <- database[[column]]
    tolerance <- 3 * getFwhm(model, masses) / (2 * sqrt(2 * log(2)))
    eligible <- eligible | vapply(seq_along(masses), function(i) {
        is.finite(masses[i]) && is.finite(tolerance[i]) &&
            any(abs(annotations$mz - masses[i]) <= tolerance[i])
    }, logical(1))
}
columns <- c("Lipid ID", "Formula (pH7.3)", "Abbreviation*", "numDoubleBond", "lipidGroup",
    "chainLength", "Exact m/z of [M-H]-", "Exact m/z of [M+H]+", "Exact m/z of [M+Na]+", "Exact m/z of [M+K]+")
saveRDS(database[eligible, columns, drop = FALSE], file.path(output, "lipids.rds"))
saveRDS(database[, columns, drop = FALSE], file.path(output, "swisslipids-species.rds"))
unlink(file.path(output, "example.log"))
saveRDS(list(source = "Cropped reduced glioblastoma example; not profile spectra",
    sourceCommit = "8b61395", sourceMd5 = unname(tools::md5sum("data/processed-example-Data.RData")),
    pixels = length(keep), bounds = lapply(coordinates[keep, ], range),
    centroided = TRUE, reduced = TRUE, normalization = "none", fwhm = "default FTICR",
    annotationFdrMaximum = max(annotations$fdr), lipidRows = sum(eligible)), file.path(output, "provenance.rds"))
cat("Native Shiny example:", length(keep), "pixels;", sum(eligible), "lipid rows\n")
print(file.info(list.files(output, full.names = TRUE))[, "size", drop = FALSE])
