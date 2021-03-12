#' Shrink data
#'
#' This function is for internal use only. It tries to shrink MSI data (in terms of size) for including it within
#' a package.
#'
#' @param x: 	   Dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param mzKeep: numeric vector, m/z values to retain.
#' @param fwhmFun: a function, fwhm as a function of m/s axis. See \code{?moleculaR::estimateFwhm}.
#' @return
#' Filtered list of \code{MALDIquant::MassPeaks} objects.
#'
#' @export
#' @keywords internal
#'
shrinkData     = function(x, mzKeep, fwhmFun) {

      if(!MALDIquant::isMassPeaksList(x)) {
            stop("Error in CreateSparseMat; x must be a list of MassPeaks objects.\n")
      }

      mzKeep      = sort(mzKeep)
      x                    = lapply(x, FUN = function(i) {

            xm            = MALDIquant::mass(i)
            #xkpIdx        = which(xm %in% mzKeep)
            xkpIdx      = MALDIquant::match.closest(x = mzKeep,
                                                    table = xm,
                                                    tolerance = fwhmFun(mzKeep))
            xkpIdx      = xkpIdx[!is.na(xkpIdx)]

            #mt            = MALDIquant::metaData(i)


            MALDIquant::createMassPeaks(mass = xm[xkpIdx],
                                        intensity = MALDIquant::intensity(i)[xkpIdx],
                                        snr = MALDIquant::snr(i)[xkpIdx],
                                        metaData = list())

      })




      return(x)

}

