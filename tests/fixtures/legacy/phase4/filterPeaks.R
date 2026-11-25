#' Filter Peaks
#'
#' This function filters out peaks which occur in less than 'minFreq' percent of all pixels.
#' This is similar to 'MALDIquant::filterPeaks' but is more memory friendly as it relies on
#' sparse data representation.
#'
#' @param x: 	Dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param minFreq: numeric falling in [0,1], peaks occuring in \code{<= minFreq * length(x)} will be omitted.
#' @return
#' Filtered list of `MassPeaks` objects.
#'
#' @export
#'
#'
filterPeaks     <- function(x, minFreq = 0.01) {

        if(!MALDIquant::isMassPeaksList(x)) {
                stop("Error in filterPeaks; x must be a list of MassPeaks objects.\n")
        }

        # first create a spars matrix per tissue region
        spData   <- createSparseMat(x)


       # find which peaks (columns) occur in less than 1% of the tissue
       relFreq              <- diff(spData$spmat@p) / nrow(spData$spmat)
       kpIdx                <- which(relFreq >= minFreq)
       .uniqueMass          <- spData$mzAxis
       kpMass               <- .uniqueMass[kpIdx]

       x                    <- lapply(x, FUN = function(i) {

              xm            <- MALDIquant::mass(i)
              xkpIdx        <- which(xm %in% kpMass)
              mt            <- MALDIquant::metaData(i)

              if(!is.null(mt$fwhm)){
                       mt$fwhm <-  mt$fwhm[xkpIdx]
              }



              MALDIquant::createMassPeaks(mass = xm[xkpIdx],
                                          intensity = MALDIquant::intensity(i)[xkpIdx],
                                          snr = MALDIquant::snr(i)[xkpIdx],
                                          metaData = mt)

       })




       return(x)

}

