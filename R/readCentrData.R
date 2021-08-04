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

      #// test if imzML of type processed, if not through an error
      # read the imzml part as a text -  this makes it faster
      txt         = readLines(con = path)

      if(any(grepl(pattern = "continuous", x = txt[3:length(txt)]))){
            stop("The input file must be of type 'processed' i.e. a centroided imzML.\n")
      }


      MALDIquantForeign::importImzMl(path = path, centroided = TRUE, verbose = verbose)

}
