.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
results <- testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)
counts <- as.data.frame(results)
write.csv(counts[, c("file", "test", "nb", "failed", "skipped", "error", "warning", "passed")],
          "docs/Cardinal-migration/test-results.csv", row.names = FALSE)
