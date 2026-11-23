# Run from the repository root with Rscript --vanilla.
lib <- file.path(getwd(), "renv", "library", "cardinal-migration")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
cat(R.version.string, "\n")
cat("Libraries:\n", paste(.libPaths(), collapse = "\n"), "\n")
packages <- c("Cardinal", "matter", "BiocParallel", "CardinalIO", "BiocManager",
              "MALDIquant", "MALDIquantForeign", "moleculaR", "testthat",
              "roxygen2", "devtools", "renv")
for (pkg in packages) {
    cat(pkg, if (requireNamespace(pkg, quietly = TRUE))
        as.character(packageVersion(pkg)) else "MISSING", "\n")
}
if (requireNamespace("Cardinal", quietly = TRUE)) {
    for (name in c("readMSIData", "MSImagingArrays", "MSImagingExperiment",
                   "mz", "intensity", "coord", "run", "spectrapply", "process",
                   "processingData", "centroided")) {
        cat("\n", name, "\n")
        if (!name %in% getNamespaceExports("Cardinal")) {
            cat("Not exported in this release\n")
            next
        }
        print(args(getExportedValue("Cardinal", name)))
        if (methods::isGeneric(name, where = asNamespace("Cardinal")))
            methods::showMethods(name, where = asNamespace("Cardinal"))
    }
}
print(sessionInfo())
