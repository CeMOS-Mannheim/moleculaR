#' Runs one of the campanion shiny apps
#'
#' Calls one of the shiny apps accompanying this package. \code{package-app} is the
#' workhorse app of the \code{moleculaR}. \code{web-app} only load an filtered
#' example dataset and is used only for showcasing the functionalities.
#'
#' @param app:    a character, name of the shiny app, \code{c("package-app", "web-app")}.
#'
#' @return
#' Has no return. To quit simply interrupt R by pressing Ctrl+c or Esc.
#'
#' @export
#'
runShiny = function(app) {
      # locate all the shiny apps that exist
      validApps = gsub(list.files(system.file("shiny", package = "moleculaR")),
                       pattern = ".R$", replacement = "")

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
      appDir = system.file("shiny", paste0(app, ".R"), package = "moleculaR")
      shiny::runApp(appDir, display.mode = "normal")
}
