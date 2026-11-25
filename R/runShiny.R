#' Run a companion Shiny application
#'
#' Launches the Cardinal-based upload/example application or the example alone.
#' MSI inputs must be a matching, known-centroided imzML/ibd pair. Spectra remain
#' file-backed in private session storage. Choose an approximate instrument
#' FWHM model or upload a separate profile spectrum TSV for measured widths;
#' reduced centroided MSI is unsuitable for peak picking or width estimation.
#' Each session supports run/ROI selection, optional Cardinal normalization,
#' single-analyte and collective maps, grouped lipid searches and portable
#' CSV/RDS/PNG downloads. Copied uploads and derived files are removed on exit.
#' The bundled example is a reduced 32-by-32 coordinate crop.
#'
#' @param app Character scalar: \code{"package-app"} (upload and example tabs)
#'   or \code{"example-app"} (example only).
#' @return Runs the application until interrupted. Press Ctrl+C or Esc to quit.
#' @export
runShiny <- function(app) {
    validApps <- c("package-app", "example-app")
    if (missing(app) || !is.character(app) || length(app) != 1L || is.na(app) || !app %in% validApps) {
        stop("Supply a valid app: 'package-app' or 'example-app'.", call. = FALSE)
    }
    appFile <- system.file("shiny", paste0(app, ".R"), package = "moleculaR", mustWork = TRUE)
    shiny::runApp(appFile, display.mode = "normal")
}
