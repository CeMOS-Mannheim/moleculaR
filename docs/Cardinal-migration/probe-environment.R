# Inspect the actual R runtime and installed MSI APIs without changing packages.
# Run from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/probe-environment.R
#
# validate-environment.R captures this output in environment.txt. Keep these
# runtime signatures distinct from the supplied older Cardinal manual.

lib <- file.path(getwd(), "renv", "library", "cardinal-migration")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
cat(R.version.string, "\n")
cat("Libraries:\n", paste(.libPaths(), collapse = "\n"), "\n")

packages <- c("Cardinal", "matter", "BiocParallel", "CardinalIO", "BiocManager", "MALDIquant",
    "MALDIquantForeign", "moleculaR", "testthat", "roxygen2", "devtools", "renv")
packages <- c(packages, "spatstat.geom", "spatstat.explore", "spatstat.random", "spatstat.utils")

for (pkg in packages) {
    cat(pkg, if (requireNamespace(pkg, quietly = TRUE))
        as.character(packageVersion(pkg)) else "MISSING", "\n")
}

# Some generic accessors are exported by their defining package rather than
# Cardinal itself. Report that distinction instead of assuming re-exports.
if (requireNamespace("Cardinal", quietly = TRUE)) {
    for (name in c("readMSIData", "MSImagingArrays", "MSImagingExperiment", "mz", "intensity",
        "coord", "run", "spectrapply", "process", "processingData", "centroided")) {
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
