#' Read Centroided imzML
#'
#' A wrapper for \code{MALDIquantForeign::importImzMl} for reading centroided (i.e. processed)
#' imzML data.
#'
#' @param path: 	character string, path to the imzML file which should be read in.
#' @param verbose:      logical, verbose output?
#'
#' @return
#' A list of \code{MALDIquant::MassPeaks} objects.
#'
#' @export
#'
#'

readCentrData     = function(path, verbose = FALSE){

      # -- To do -- #

      #// test if imzML of type processed, if not through an error, use XML package

      # ---------- #


      MALDIquantForeign::importImzMl(path = path, centroided = TRUE, verbose = verbose)

}
