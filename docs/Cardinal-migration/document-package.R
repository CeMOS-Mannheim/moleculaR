.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
# Regenerate the migration API and namespace without unrelated Rd formatting churn.
document_migration <- function() {
    unrelated <- setdiff(list.files("man", full.names = TRUE),
                         file.path("man", c("readMsi.Rd", "searchAnalyte.Rd")))
    original <- lapply(unrelated, function(path) readBin(path, "raw", file.size(path)))
    on.exit(for (i in seq_along(unrelated)) writeBin(original[[i]], unrelated[i]))
    roxygen2::roxygenize(".", roclets = c("rd", "namespace", "collate"))
}
document_migration()
