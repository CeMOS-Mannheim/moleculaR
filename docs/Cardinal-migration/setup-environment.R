# Isolated library: never restore the old project lockfile or update user libraries.
# Run from the repository root with Rscript --vanilla.
lib <- file.path(getwd(), "renv", "library", "cardinal-migration")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 300)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = lib, type = "binary")
options(repos = BiocManager::repositories(version = "3.23"))
packages <- c("BiocVersion", "Cardinal", "CardinalIO", "matter", "BiocParallel",
              "MALDIquant", "MALDIquantForeign", "import", "spatstat", "Matrix",
              "raster", "viridis", "vioplot", "shiny", "shinyWidgets",
              "shinythemes", "testthat", "roxygen2", "pkgload", "renv")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, type = "binary")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
source("docs/Cardinal-migration/probe-environment.R")
