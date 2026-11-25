#' Filter Peaks
#'
#' This function filters out peaks which occur in less than 'minFreq' percent of all pixels.
#' This is similar to 'MALDIquant::filterPeaks' but is more memory friendly as it relies on
#' sparse data representation.
#'
#' @param x An aligned Cardinal MSImagingExperiment with known centroided status,
#'   or a temporary legacy MassPeaks list. Native arrays require explicit alignment.
#' @param minFreq Numeric between zero and one. Keep features with frequency >= minFreq,
#'   including equality.
#' @param featureChunkSize Native input only: optional number of feature rows
#'   per block. Default targets Cardinal's configured byte budget, or 16 MiB.
#' @details Native frequencies count positive intensities per acquired pixel.
#'   Stored zeros are absent; duplicate shared rows remain separate features.
#'   Counts and frequencies are recomputed in featureData. Coordinates, metadata,
#'   empty feature results and completed processing history are preserved.
#'   Complete existing processing queues with Cardinal::process() first.
#' @return
#' A filtered native Cardinal object, or the legacy filtered MassPeaks list.
#'
#' @export
#' @md
#'
#'
filterPeaks     <- function(x, minFreq = 0.01, featureChunkSize = NULL) {

        if (.msiIsObject(x)) return(.filterCardinalPeaks(x, minFreq, featureChunkSize))
        if (!is.null(featureChunkSize)) stop("featureChunkSize requires Cardinal input.")

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

