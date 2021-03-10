#' Normalize Intensity
#'
#' This normalizes spectral intensities such that for every spectrum a median
#' fold-change of peak intensities is calculated with reference to the mean spectrum
#' and normalizes based on that.
#'
#' @param x: 	        a list of \code{MassPeaks} objects.
#' @param meanMethod:   see \code{?MALDIquant::mergeMassPeaks} for more info.
#' @return
#' Intensity-normalized list of \code{MassPeaks} objects.
#'
#' @references
#' Veselkov, Kirill, et al. "BASIS: High-performance bioinformatics platform for processing of large-scale mass
#' spectrometry imaging data in chemically augmented histology." Scientific reports 8.1 (2018): 1-11.
#'
#' @export
#'
foldChangeNorm          = function(x, meanMethod = "mean") {


        # first find the mean mass Peaks object
        dataMean        = MALDIquant::mergeMassPeaks(x, method = meanMethod)


        # for every spectrum in each tissue section find median fold-change to the mean spectrum and normalize

        x               = lapply(x, FUN = function(spect) {


                mSpect  = MALDIquant::mass(spect)
                mMean   = MALDIquant::mass(dataMean)
                iSpect  = MALDIquant::intensity(spect)
                iMean   = MALDIquant::intensity(dataMean)

                idxInMean = which(mMean %in% mSpect)
                idxInSpect= which(mSpect %in% mMean[idxInMean])

                fc      = iSpect[idxInSpect] / iMean[idxInMean]
                medfc   = median(fc)

                iSpect  = iSpect / medfc

                MALDIquant::intensity(spect) = iSpect



                return(spect)

       })

       return(x)




}

