.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
options(repos = BiocManager::repositories(version = "3.23"), timeout = 300)
report <- capture.output({
    source("docs/Cardinal-migration/probe-environment.R")
    cat("\nBioconductor version:", as.character(BiocManager::version()), "\n")
    validity <- BiocManager::valid(lib.loc = .libPaths(), checkBuilt = TRUE)
    print(validity)
})
writeLines(sub("[ \t]+$", "", report), "docs/Cardinal-migration/environment.txt", useBytes = TRUE)
if (!isTRUE(validity)) stop("Resolve Bioconductor validity issues before snapshotting.")
Sys.setenv(RENV_CONFIG_CACHE_ENABLED = "FALSE",
           RENV_PATHS_ROOT = file.path(getwd(), "renv", "local", "migration-state"))
renv::snapshot(project = getwd(), library = .libPaths(), type = "all", prompt = FALSE)
# Official Windows binaries may carry an r-universe Repository metadata field.
# Pin their source to the selected official Bioconductor release for restore.
bioc <- c("Biobase", "BiocGenerics", "BiocParallel", "BiocVersion", "Cardinal",
          "CardinalIO", "ProtGenerics", "S4Vectors", "matter")
records <- lapply(bioc, function(pkg) list(Package = pkg,
    Version = as.character(packageVersion(pkg)), Source = "Bioconductor",
    Repository = "Bioconductor"))
names(records) <- bioc
renv::record(records, project = getwd())
