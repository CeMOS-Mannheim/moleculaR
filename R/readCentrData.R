#' Read imzML data
#'
#' A wrapper for \code{MALDIquantForeign::importImzMl} for reading centroided (i.e. processed)
#' and continuous imzML data.
#'
#' @param path: 	character string, path to the imzML file which should be read in.
#' @param verbose:      logical, verbose output?
#' @param ...: other arguments passed to \code{MALDIquantForeign::importImzMl}.
#'
#' @return
#' A list of \code{MALDIquant::MassPeaks} objects if the input file is of 'processed' type or
#' a list of \code{MALDIquant::MassSpectrum} objects if the input file is of 'continuous' type.
#'
#' @export
#'
#'

readCentrData     = function(path, verbose = FALSE, ...){

      #// test if imzML of type processed, if not through an error
      # read the imzml part as a text -  this makes it faster
      txt         = readLines(con = path, n = 50)

      if(any(grepl(pattern = "continuous", x = txt[3:length(txt)]))){

         cat("Input dataset is of type 'continuous'. This could lead to high memory consumption. \n")

         MALDIquantForeign::importImzMl(path = path, centroided = FALSE, verbose = verbose, ...)

      } else {

         cat("Input dataset is of type 'processed'.\n")

         MALDIquantForeign::importImzMl(path = path, centroided = TRUE, verbose = verbose, ...)

      }


}
