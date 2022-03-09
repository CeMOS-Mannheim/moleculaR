#' Normalize Intensity
#'
#' This normalizes spectral intensities based on a number of different methods. such that for every spectrum a median
#' fold-change of peak intensities is calculated with reference to the mean spectrum
#' which is then taken as a normalization factor for that particular spectrum.
#'
#' @param x:               a list of \code{MassPeaks} objects.
#' @param method:          the normalization method,  has to be one of c('medianFoldChange', 'RMS', 'TIC', 'median').
#' @param meanMethod:   only used for `medianFoldChange` normalization. See
#' \code{?MALDIquant::mergeMassPeaks} for more info.
#'
#' @details
#' The `medianFoldChange` normalization is done such that for every spectrum a median
#' fold-change of peak intensities is calculated with reference to the mean spectrum
#' which is then taken as a normalization factor for that particular spectrum (see
#' reference [1]). The `RMS`, `TIC` and `median` normalization methods are
#' equivalent to their counterparts in other packages (ex. MALDIquant), they are
#' just adapted to process centroided data of `MassPeak` Objects. Note that for
#' `medianFoldChange` normalization, it is important to run peak binning beforehand.
#'
#' @return
#' Intensity-normalized list of \code{MassPeaks} objects.
#'
#' @references
#' 1. Veselkov, Kirill, et al. "BASIS: High-performance bioinformatics platform for processing of large-scale mass
#'    spectrometry imaging data in chemically augmented histology." Scientific reports 8.1 (2018): 1-11.
#' 2. Deininger, Sören-Oliver, et al. "Normalization in MALDI-TOF imaging datasets of proteins: practical
#'    considerations." Analytical and bioanalytical chemistry 401.1 (2011): 167-181.
#'
#'
#' @export
#'

normalizeIntensity <- function(x, method, meanMethod = "mean"){


         x <- switch(method,
                     "medianFoldChange" = .foldChangeNorm(x, meanMethod),
                     "RMS" = .rmsNorm(x),
                     "TIC" = .ticNorm(x),
                     "median" = .medianNorm(x),
                     stop("method has to be one of ",
                          "c('medianFoldChange', 'RMS', 'TIC', 'median').\n"))

         return(x)

}


#' Normalize Intensity based on foldchange
#'
#' This normalizes spectral intensities such that for every spectrum a median
#' fold-change of peak intensities is calculated with reference to the mean spectrum
#' which is then taken as a normalization factor for that particular spectrum.
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
#' @keywords internal
#'

.foldChangeNorm          <- function(x, meanMethod = "mean") {


         # first find the mean mass Peaks object
         dataMean        <- MALDIquant::mergeMassPeaks(x, method = meanMethod)


         # for every spectrum in each tissue section find median fold-change to the mean spectrum and normalize

         x               <- lapply(x, FUN = function(spect) {


                  mSpect  <- MALDIquant::mass(spect)
                  mMean   <- MALDIquant::mass(dataMean)
                  iSpect  <- MALDIquant::intensity(spect)
                  iMean   <- MALDIquant::intensity(dataMean)

                  idxInMean <- which(mMean %in% mSpect)
                  idxInSpect<- which(mSpect %in% mMean[idxInMean])

                  fc      <- iSpect[idxInSpect] / iMean[idxInMean]
                  medfc   <- median(fc)

                  iSpect  <- iSpect / medfc

                  MALDIquant::intensity(spect) <- iSpect



                  return(spect)

         })

         return(x)




}




#' RMS Intensity Normalization
#'
#' This normalizes spectral intensities to their RMS.
#'
#' @param x: 	      a list of `MassPeaks` or `MassSpectrum` objects.
#' @return
#' RMS intensity-normalized list of `MassPeaks` or `MassSpectrum` objects.
#'
#'  @keywords internal
#'

.rmsNorm            <- function(x) {


         x           <- lapply(x, FUN = function(spect) {


                  i  <- MALDIquant::intensity(spect)

                  MALDIquant::intensity(spect) <- i / sqrt(sum(i^2)/(length(i)))

                  return(spect)

         })

         return(x)


}


#' TIC Intensity Normalization
#'
#' This normalizes spectral intensities to their TIC.
#'
#' @param x: 	      a list of `MassPeaks` or `MassSpectrum` objects.
#' @return
#' RMS intensity-normalized list of `MassPeaks` or `MassSpectrum` objects.
#'
#' @keywords internal
#'

.ticNorm            <- function(x) {


         x           <- lapply(x, FUN = function(spect) {


                  i  <- MALDIquant::intensity(spect)

                  MALDIquant::intensity(spect) <- i / sum(i)

                  return(spect)

         })

         return(x)


}


#' Median Intensity Normalization
#'
#' This normalizes spectral intensities to their median.
#'
#' @param x: 	      a list of `MassPeaks` or `MassSpectrum` objects.
#' @return
#' RMS intensity-normalized list of `MassPeaks` or `MassSpectrum` objects.
#'
#'  @keywords internal
#'

.medianNorm            <- function(x) {


         x           <- lapply(x, FUN = function(spect) {


                  i  <- MALDIquant::intensity(spect)

                  MALDIquant::intensity(spect) <- i / median(i)

                  return(spect)

         })

         return(x)


}

