#' Runs one of the campanion shiny apps
#'
#' Calls one of the shiny apps accompanying this package.
#'
#' @param app:    a character, name of the shiny app, c("package-app", "web-app").
#'
#' @return
#' A sparse matrix, see \code{?Matrix::sparseMatrix} for more info.
#'
#' @export
#'
runShiny = function(app) {
      # locate all the shiny apps that exist
      validApps = list.files(system.file("shiny", package = "moleculaR"))

      validAppsMsg =
            paste0(
                  "Valid shiny apps are: '",
                  paste(validApps, collapse = "', '"),
                  "'")

      # if an invalid app is given, throw an error
      if (missing(app) || !nzchar(app) ||
          !(app %in% validApps)) {
            stop(
                  'Please run `runShiny()` with a valid app as an argument.\n',
                  validAppsMsg,
                  call. = FALSE)
      }

      # find and launch the app
      appDir = system.file("shiny-examples", app, package = "moleculaR")
      shiny::runApp(appDir, display.mode = "normal")
}
