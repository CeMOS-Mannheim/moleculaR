# Regenerate API help, exports, and source collation with roxygen2.
# Usage from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/document-package.R
#
# Preserve unrelated Rd files byte-for-byte. An upgrade of roxygen2 should not
# turn a focused migration change into a package-wide documentation reformat.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
document_migration <- function() {
    unrelated <- setdiff(list.files("man", full.names = TRUE), file.path("man", c("readMsi.Rd",
        "searchAnalyte.Rd", "batchAnalyteSearch.Rd", "batchLipidSearch.Rd",
        "readSingleSpect.Rd", "estimateFwhm.Rd", "useDefaultFwhm.Rd",
        "normalizeIntensity.Rd", "filterPeaks.Rd")))
    original <- lapply(unrelated, function(path) readBin(path, "raw", file.size(path)))
    on.exit(for (i in seq_along(unrelated)) writeBin(original[[i]], unrelated[i]))

    roxygen2::roxygenize(".", roclets = c("rd", "namespace", "collate"))
}

document_migration()
