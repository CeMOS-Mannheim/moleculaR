#' Normalize Intensity
#'
#' This normalizes spectral intensities based on a number of different methods.
#'
#' @param x A native Cardinal object or a temporary legacy MassPeaks list.
#' @param method The normalization method, one of
#' c('medianFoldChange', 'RMS', 'TIC', 'median', 'IS').
#' @param meanMethod Reference aggregation for medianFoldChange: mean, median or sum.
#' @param mz A numeric, the m/z value of the internal standard (IS) -
#' only relevant for `method='IS'`.
#' @param tol A numeric, the absolute tolerance for m/z (in Da) -
#' only relevant for `method='IS'`.
#' @param scale Native input only: positive output scale, or NA for Cardinal's
#'   defaults (maximum spectrum length/shared feature count for TIC/RMS, one
#'   otherwise). Use one for division-only normalization.
#' @param zeroAsMissing Native medianFoldChange only: interpret shared zeros as
#'   absent centroid peaks (TRUE), or explicitly observed zeros (FALSE).
#' @param outfile Native input only: new intensity backing file. Defaults to
#'   temporary session storage. Keep this file and the original m/z source files
#'   accessible, or provide a persistent destination. Existing paths are refused.
#' @param featureChunkSize Optional feature rows per reference block for native
#'   medianFoldChange. The default targets Cardinal's byte budget, or 16 MiB.
#' @param BPPARAM Optional native BiocParallel backend; serial by default.
#' @param verbose Logical; show native processing progress.
#'
#' @details
#' Native TIC, RMS and reference scaling use the same exported matter rescaling
#' functions as Cardinal. Empty/all-zero spectra remain zero. Non-finite inputs
#' and undefined nonzero-spectrum factors raise errors. IS uses the nearest
#' reference within tol on a sorted mass axis; absent/zero standards yield zero
#' spectra, following Cardinal's zero-reference policy. It does not apply the
#' legacy median imputation. Default TIC/RMS scales differ from legacy values.
#' Factors are recorded in pixelData and completed processing in Cardinal history.
#' Existing queued preprocessing must be processed explicitly before this call.
#'
#' Native medianFoldChange requires aligned centroided MSImagingExperiment data.
#' Reference-zero features are excluded from ratios; zeroAsMissing controls the
#' interpretation of other zeros. Arrays are never aligned implicitly.
#'
#' On the temporary legacy path, `medianFoldChange` is done such that for every spectrum a median
#' fold-change of peak intensities is calculated with reference to the mean spectrum
#' which is then taken as a normalization factor for that particular spectrum (see
#' reference 1). The `RMS`, `TIC` and `median` normalization methods are
#' equivalent to their counterparts in other packages (ex. MALDIquant), they are
#' just adapted to process centroided data of `MassPeak` Objects. Note that for
#' `medianFoldChange` normalization, it is important to run peak binning beforehand.
#' The `IS` method assumes the presence of an internal standard at `mz` m/z against
#' which all intensities are normalized (pixel-wise).
#'
#' @return
#' Native input returns the same Cardinal class with completed, file-backed
#' intensities in a new backing file. Legacy input returns a normalized list.
#'
#' @references
#' 1. Veselkov, Kirill, et al. "BASIS: High-performance bioinformatics platform for processing of large-scale mass
#'    spectrometry imaging data in chemically augmented histology." Scientific reports 8.1 (2018): 1-11.
#' 2. Deininger, Sören-Oliver, et al. "Normalization in MALDI-TOF imaging datasets of proteins: practical
#'    considerations." Analytical and bioanalytical chemistry 401.1 (2011): 167-181.
#'
#'
#' @export
#' @md
#'

normalizeIntensity <- function(x, method, meanMethod = "mean", mz = NULL, tol = NULL,
                               scale = NA_real_, zeroAsMissing = TRUE, outfile = NULL,
                               featureChunkSize = NULL, BPPARAM = NULL, verbose = FALSE) {

         if (.msiIsObject(x)) {
                  return(.normalizeCardinal(x, method, meanMethod, mz, tol, scale,
                         zeroAsMissing, outfile, featureChunkSize, BPPARAM, verbose))
         }
         if (!missing(scale) || !missing(zeroAsMissing) || !is.null(outfile) ||
             !is.null(featureChunkSize) || !is.null(BPPARAM) || isTRUE(verbose)) {
                  stop("Native normalization options require Cardinal input.")
         }


         if(method == "IS"){
                  if(is.null(mz))
                           stop("'mz' must be provided for 'method=IS'.\n")

                  if(is.null(tol))
                           stop("'tol' must be provided for 'method=IS'.\n")
         }

         x <- switch(method,
                     "medianFoldChange" = .foldChangeNorm(x, meanMethod),
                     "RMS" = .rmsNorm(x),
                     "TIC" = .ticNorm(x),
                     "median" = .medianNorm(x),
                     "IS" = .isNorm(x, mz, tol),
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


#' IS Intensity Normalization
#'
#' This normalizes spectral intensities to a provided M/z value i.e. an intenanl
#' standard (IS).
#'
#' @param x: 	      a list of `MassPeaks` or `MassSpectrum` objects.
#' @param mz:           a numeric, the m/z value of the IS.
#' @param tol:          a numeric, the absolute tolerance for m/z (in Da).
#' @return
#' IS intensity-normalized list of `MassPeaks` or `MassSpectrum` objects.
#'
#' @keywords internal
#'

.isNorm            <- function(x, mz, tol) {


         # find the intensity of IS in each pixel
         is <- .findis(x, mz, tol)

         # replace NAs(no detection of IS) with the median of all IS values
         is[is.na(is)] <- median(is, na.rm = TRUE)

         x           <- lapply(seq_along(x), FUN = function(i) {

                  int  <- MALDIquant::intensity(x[[i]])

                  MALDIquant::intensity(x[[i]]) <- int / is[i]

                  return(x[[i]])

         })

         return(x)


}

.findis <- function(x, mz, tol){

         is <- sapply(x, FUN = function(spect) {

                  idx <- MALDIquant::match.closest(x = mz,
                                                   table = spect@mass,
                                                   tolerance = tol)

                  if(is.na(idx)){
                           return(idx)
                  } else {
                           return(spect@intensity[idx])
                  }

         })

         return(is)

}
