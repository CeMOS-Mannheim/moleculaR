#' Read an MSI dataset using Cardinal
#'
#' Import an imzML/ibd pair as a native Cardinal object. Spectral data remain
#' file-backed by default. Keep the original files accessible for the lifetime
#' of the object, including after saving or restoring it with R serialization.
#'
#' @param path Path to an imzML file with a paired ibd file of the same stem.
#' @param memory Logical; load spectra into memory if TRUE (default FALSE).
#' @param as Cardinal representation: \code{"MSImagingArrays"} (default),
#'   \code{"MSImagingExperiment"}, or \code{"auto"}. Requesting a shared-axis
#'   representation for unaligned spectra may align or resample them.
#' @param run Optional single run name to retain. Other runs remain separate
#'   spatial domains. If NULL, import all runs.
#' @param ... Arguments to \code{Cardinal::readMSIData()}, including
#'   \code{check = TRUE} for binary UUID/checksum validation. With this option,
#'   Cardinal reader warnings are treated as errors.
#'
#' @details Import preserves spectrum representation and metadata without
#'   normalization, peak picking, or mandatory alignment. Profile data and unknown
#'   centroid status may be imported; \code{searchAnalyte()} requires known
#'   centroided spectra and completed preprocessing. Select a single run before
#'   constructing a probability map. For an existing arrays object use
#'   \code{x[Cardinal::run(x) == "sample"]}; for an experiment use
#'   \code{x[, Cardinal::run(x) == "sample", drop = FALSE]}.
#'   The supported migration environment is R 4.6 / Bioconductor 3.23 with
#'   Cardinal 3.14.0. See the migration documentation for reproducible setup.
#' @return A native \code{MSImagingArrays} or \code{MSImagingExperiment} object.
#' @seealso \code{\link{searchAnalyte}}, \code{Cardinal::readMSIData}
#' @export
readMsi <- function(path, memory = FALSE, as = "MSImagingArrays", run = NULL, ...) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
        stop("path must be a single imzML filename.")
    if (!file.exists(path) || isTRUE(file.info(path)$isdir))
        stop("The imzML file does not exist: ", path)
    if (tolower(tools::file_ext(path)) != "imzml")
        stop("readMsi() currently supports imzML files only.")
    ibd <- paste0(tools::file_path_sans_ext(path), ".ibd")
    if (!file.exists(ibd) || isTRUE(file.info(ibd)$isdir))
        stop("The paired ibd file does not exist: ", ibd)
    if (!is.logical(memory) || length(memory) != 1L || is.na(memory))
        stop("memory must be TRUE or FALSE.")
    as <- match.arg(as, c("MSImagingArrays", "MSImagingExperiment", "auto"))
    dots <- list(...)
    if (isTRUE(dots$parse.only))
        stop("readMsi() returns MSI objects; use Cardinal::readMSIData() for parse.only.")
    x <- withCallingHandlers(
        Cardinal::readMSIData(file = path, memory = memory, as = as, ...),
        warning = function(w) {
            if (isTRUE(dots$check)) stop(conditionMessage(w), call. = FALSE)
        })
    if (!is.null(run)) {
        if (!is.character(run) || length(run) != 1L || is.na(run) || !nzchar(run))
            stop("run must be a single run name.")
        pixels <- .msiRun(x) == run
        if (!any(pixels)) stop("Requested run was not found: ", run)
        x <- .msiSelectPixels(x, pixels)
    }
    .msiValidate(x)
    x
}
