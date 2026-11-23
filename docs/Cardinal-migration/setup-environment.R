# Install development dependencies into a repository-local library.
# Run from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/setup-environment.R
#
# This convenience installer uses available binaries for Bioconductor 3.23.
# Use renv.lock to restore the exact versions recorded by validation. No package
# is installed into or updated inside the user's global library.

# Keep the library inside ignored renv storage, separate from the old runtime.
lib <- file.path(getwd(), "renv", "library", "cardinal-migration")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

# Obtain the Bioconductor repository resolver before selecting that release.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 300)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib,
    type = "binary")
options(repos = BiocManager::repositories(version = "3.23"))

# List direct runtime dependencies and the tools used by the validation scripts.
# Spatial functions use the child packages directly; the former import and
# spatstat meta-package dependencies are no longer needed.
packages <- c(
    # Bioconductor MSI infrastructure.
    "BiocVersion", "Cardinal", "CardinalIO", "matter", "BiocParallel",

    # Existing spectral reference and spatial analysis.
    "MALDIquant", "MALDIquantForeign", "Matrix",
    "spatstat.geom", "spatstat.explore", "spatstat.random", "spatstat.utils",

    # Existing visualization/application dependencies.
    "raster", "viridis", "vioplot", "shiny", "shinyWidgets", "shinythemes",

    # Development and reproducibility tools.
    "testthat", "roxygen2", "pkgload", "renv"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
    install.packages(missing, lib = lib, type = "binary")
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))

# Print the versions and actual Cardinal signatures after installation.
source("docs/Cardinal-migration/probe-environment.R")
