# Build and check the package in isolated subprocesses using the local library.
# Usage from the repository root:
#   Rscript --vanilla docs/Cardinal-migration/check-package.R
#
# Build artifacts stay in ignored renv/local storage. Compact UTF-8 reports are
# written beside this script; existing manual regression scripts are exercised
# by R CMD check as well as the testthat suite. PDF manuals and vignette rebuilds
# are disabled, so this does not replace the final release validation protocol.

.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
root <- normalizePath(getwd(), winslash = "/")
output <- file.path(root, "renv/local/cardinal-check")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Normalize Windows subprocess line endings for readable, diffable reports.
write_report <- function(text, path) {
    lines <- unlist(strsplit(gsub("\r", "", text, fixed = TRUE), "\n", fixed = TRUE), use.names = FALSE)
    lines <- sub("[ \t]+$", "", lines)
    while (length(lines) && tail(lines, 1L) == "") {
        lines <- head(lines, -1L)
    }
    writeLines(lines, path)
}

# The host's C.UTF-8 locale is not recognized by Windows R. Use C only inside
# these subprocesses, without changing the user's system or session settings.
checkenv <- c(
    RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
    LC_ALL = "C",
    LANG = "C",
    LC_CTYPE = "C",
    LC_COLLATE = "C",
    LC_MONETARY = "C",
    LC_TIME = "C",
    `_R_CHECK_FORCE_SUGGESTS_` = "false"
)

# Check the built source archive, not an accidentally installed older package.
build <- callr::rcmd("build", c("--no-build-vignettes", "--no-manual", root), wd = output,
    libpath = .libPaths(), fail_on_status = FALSE, env = checkenv)
write_report(c(build$stdout, build$stderr), "docs/Cardinal-migration/package-build.txt")
cat(build$stdout, build$stderr, sep = "\n")
if (build$status != 0L) stop("Package build failed.")
archive <- list.files(output, "moleculaR_.*[.]tar[.]gz$", full.names = TRUE)
archive <- archive[which.max(file.info(archive)$mtime)]

check <- callr::rcmd("check", c("--no-manual", "--no-build-vignettes", archive), wd = output,
    libpath = .libPaths(), fail_on_status = FALSE, env = checkenv)
write_report(c(check$stdout, check$stderr), "docs/Cardinal-migration/package-check.txt")

cat(tail(strsplit(check$stdout, "\n", fixed = TRUE)[[1]], 20), sep = "\n")
quit(status = check$status)
