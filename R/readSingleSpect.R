#' Read a single spectrum
#'
#' This reads a single spectrum from a tsv file containing two columnd for m/z and intensity,
#' respectively.
#'
#' @param path: 	character string, path to the imzML file which should be read in.
#'
#' @return
#' A \code{MALDIquant::MassSpectrum} object.
#'
#' @export
#'
#'

readSingleSpect     = function(path){

      #// read tsv
      df          = read.delim(file = path)

      #// convert to Massspectrum class
      ms          = MALDIquant::createMassSpectrum(mass = df[[1]], intensity = df[[2]])

      return(ms)

}
