.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
root <- normalizePath(getwd(), winslash = "/")
output <- file.path(root, "renv/local/cardinal-check")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_report <- function(text, path) {
    lines <- unlist(strsplit(gsub("\r", "", text, fixed = TRUE), "\n", fixed = TRUE),
                    use.names = FALSE)
    lines <- sub("[ \t]+$", "", lines)
    while (length(lines) && tail(lines, 1L) == "") lines <- head(lines, -1L)
    writeLines(lines, path)
}
checkenv <- c(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE", LC_ALL = "C", LANG = "C",
              LC_CTYPE = "C", LC_COLLATE = "C", LC_MONETARY = "C", LC_TIME = "C",
              "_R_CHECK_FORCE_SUGGESTS_" = "false")
build <- callr::rcmd("build", c("--no-build-vignettes", "--no-manual", root),
    wd = output, libpath = .libPaths(), fail_on_status = FALSE,
    env = checkenv)
write_report(c(build$stdout, build$stderr), "docs/Cardinal-migration/package-build.txt")
cat(build$stdout, build$stderr, sep = "\n")
if (build$status != 0L) stop("Package build failed.")
archive <- list.files(output, "moleculaR_.*[.]tar[.]gz$", full.names = TRUE)
archive <- archive[which.max(file.info(archive)$mtime)]
check <- callr::rcmd("check", c("--no-manual", "--no-build-vignettes", archive),
    wd = output, libpath = .libPaths(), fail_on_status = FALSE,
    env = checkenv)
write_report(c(check$stdout, check$stderr), "docs/Cardinal-migration/package-check.txt")
cat(tail(strsplit(check$stdout, "\n", fixed = TRUE)[[1]], 20), sep = "\n")
quit(status = check$status)
