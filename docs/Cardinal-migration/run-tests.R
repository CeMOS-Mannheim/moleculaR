# Run the automatic regression suite against the current source checkout.
# Usage from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/run-tests.R
#
# Tests stop this script on failure. The CSV records counts per test, including
# warnings/skips, so the change log can cite concrete validation results.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)

results <- testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)

counts <- as.data.frame(results)
write.csv(counts[, c("file", "test", "nb", "failed", "skipped", "error", "warning", "passed")],
    "docs/Cardinal-migration/test-results.csv", row.names = FALSE)
