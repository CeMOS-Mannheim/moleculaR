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
#' A list of \code{MALDIquant::MassPeaks} objects if the input file is of 'processed' type.
#'
#' @export
#'
#'

readCentrData     = function(path, verbose = FALSE, ...){

      #// test if imzML of type processed, if not through an error
      # read the imzml part as a text -  this makes it faster
      txt         = readLines(con = path, n = 50)

      if(any(grepl(pattern = "continuous", x = txt[3:length(txt)]))){

         stop("Input dataset is of type 'continuous'. Only 'processed' type is supported so far. \n")

      } else {

         cat("Input dataset is of type 'processed'.\n")

         MALDIquantForeign::importImzMl(path = path, centroided = TRUE, verbose = verbose, ...)

      }


}
