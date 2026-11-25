#' Read a single spectrum
#'
#' This reads a single spectrum from a tsv file containing two columnd for m/z and intensity,
#' respectively. Note that the first row is counted as the header.
#'
#' @param path Character string, path to the two-column TSV file.
#' @param header Logical; whether the first row is a header (default TRUE,
#'   preserving the historical reader). Use FALSE for headerless numeric files.
#'
#' @return
#' A numeric \code{moleculaRSpectrum} list with \code{mz} and \code{intensity}
#' vectors, accepted directly by \code{estimateFwhm()}. Masses must be positive
#' and strictly increasing; intensities must be finite and non-negative.
#'
#' @export
#'
#'

readSingleSpect <- function(path, header = TRUE) {
      if (!is.character(path) || length(path) != 1L || is.na(path) ||
          !file.exists(path) || isTRUE(file.info(path)$isdir)) stop("path must name an existing spectrum TSV file.")
      if (!is.logical(header) || length(header) != 1L || is.na(header)) stop("header must be TRUE or FALSE.")
      df <- utils::read.delim(path, header = header, row.names = NULL,
                            stringsAsFactors = FALSE, check.names = FALSE)
      if (ncol(df) != 2L) stop("The spectrum TSV must contain exactly two columns.")
      .numericSpectrum(df[[1]], df[[2]])
}
