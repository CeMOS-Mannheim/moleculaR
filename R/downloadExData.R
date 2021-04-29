#' Download example MSI data
#'
#' This downloads an example dataset for package illustration and development.
#'
#' @param path:   a character string, path to the directory where the downloaded
#' file is saved. Defaults to \code{tempdir()}. The dataset is save as
#' \code{moleculaR-example-data.imzML}.
#' @return
#' Returns the path to the downloaded example data file.
#'
#' @export
#'
#'
downloadExData    = function(path = tempdir()) {

      #// imzml + ibd
      imzmlUrl    = scan(system.file("extdata", "dataUrl-imzml.txt", package = "moleculaR", mustWork = TRUE),
                         what = "character", sep = " ")
      ibdUrl      = scan(system.file("extdata", "dataUrl-ibd.txt", package = "moleculaR", mustWork = TRUE),
                         what = "character", sep = " ")

      #// download
      download.file(url = imzmlUrl, destfile = file.path(path, "moleculaR-example-data.imzML"))
      download.file(url = ibdUrl, destfile = file.path(path, "moleculaR-example-data.ibd"))

      return(file.path(path, "moleculaR-example-data.imzML"))

}

