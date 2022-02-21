#' Estimation of fwhm as a function of m/z axis
#'
#' This creates a linear interpolator function that approximats fwhm as a function
#' of m/z for the given single spectrum or a list of mass spectra.
#'
#' @param s: a single profile (continuous) spectrum of type `MALDIquant::MassSpectrum` of a list thereof. This
#' could also be a list of `MALDIquant::MassPeaks` object but with conditions, see Details.
#' @param spectraSampling:  an integer specifying how many spectra should be used for fwhm estimation
#' (the lower the faster) when `s` is a list of `MALDIquant::MassSpectrum` objects. Defaults to 99.
#' @param peakSampling:  an integer specifying how many detected peaks (in each spectrum) to consider for fwhm estimation,
#' the lower the faster. Defaults to 99.
#' @param numCores: an integer specifying the number of cores used for the FWHM computations when s is a
#' `MassSpectrum` or a list thereof. This is ignored on windows.
#' @param plot:      whether to plot the result.
#' @param savePlot:  either \code{NULL} or file path to save a plot as svg.
#'
#' @details
#' The input `s` could be a single profile (continuous) `MassSpectrum` object. This could be supplied externally
#' or via `readSingleSpect` function (see `?readSingleSpect`). This is normally the case when dealing with MRMS (FTICR)
#' centroided data which does not contain any information about peaks FWHM values (for example when read from a centroided
#' imzML file).
#'
#' `s` could also be a list of `MassSpectrum` objects. This could be the case when dealing with continuous TOF data.
#' In this case, the list of spectra is sampled and FWHM values are computed based on these.
#'
#' If `s` is a list of centroided `MassPeaks` objects, the function assumes that FWHM values have been
#' (externally) pre-computed and attached to each centroided spectrum (`MassPeaks` object) via the
#' `metaData` slot (see `?MALDIquant::MassPeaks`) such that for every centroided spectrum `s[[i]]` the
#' FWHM values corresponding to all detected peaks `s[[i]]@mass` are stored in `s[[i]]@metaData$fwhm`.
#'
#' @return
#' Returns an S3 object 'fwhm' containing a linear interplotor function 'fwhmInterpolator' in addition to
#' m/z values of the peaks and their corresponding fwhm values.
#'
#' @export
#' @include manualSpatstatImport.R
#'
estimateFwhm            <- function(s, spectraSampling = 10, peakSampling = 1000,
                                    numCores = 1, plot = FALSE, savePlot = NULL) {


      if(is.list(s)){

            if(MALDIquant::isMassSpectrumList(s)) {


                  #// sample spectra - here s is a list of MassPSpectrum objects
                  spectraSampling <- ifelse(length(s) < spectraSampling, length(s), spectraSampling)
                  idx         <- sample(x = seq_along(s), size = spectraSampling)
                  s           <- s[idx]

                  #// detect peaks - here p is a list of MassPeaks objects
                  p           <- MALDIquant::detectPeaks(object = s, method = "SuperSmoother", SNR = 3)

                  #// override numCores if OS is windows
                  numCores    <- ifelse(.Platform$OS.type == "windows", 1, numCores)

                  #// compute fwhm at every peak in each spectrum - time consuming step
                  fwhmdf      <- parallel::mclapply(X = seq_along(s), mc.cores = numCores,
                                                    FUN = function(i){
                                                          .getFwhm(s[[i]], .samplep(p[[i]], peakSampling))
                                                    })

                  #// put everything together
                  fwhmdf      <- do.call("rbind", fwhmdf)
                  fwhmdf      <- fwhmdf[order(fwhmdf$peaks), ]

            } else{

                  # This part assumes the input is a list of MassPeaks objects
                  # with fwhm values precomputed and attached to the 'metaData'
                  # slot of each spectrum on MALDIquant::metaData(spectrum)$fwhm

                  if(!MALDIquant::isMassPeaksList(s)){
                        stop("The input s is not a list of MassSpectrum/MassPeaks objects. \n")
                  }

                  if(is.null(MALDIquant::metaData(s[[1]])$fwhm)){
                        stop("The input s is a list of MassPeaks but ",
                             "does not contain any precomputed ",
                             "fwhm data, see details in ?estimateFWHM. \n")
                  }

                  if(length(s[[1]]@mass) != length(s[[1]]@metaData$fwhm)){
                        stop("Number of peaks does not equal number of FWHM values ",
                             "attached to s. \n")
                  }

                  spectraSampling <- ifelse(length(s) < spectraSampling, length(s), spectraSampling)
                  idx         <- sample(x = seq_along(s), size = spectraSampling)
                  s           <- s[idx]

                  #// collect fwhm and corresponding masses.
                  fwhmdf       <- lapply(s, FUN = function(i) {

                        isampled <- .samplep(i, peakSampling)

                        data.frame(peaks = MALDIquant::mass(isampled),
                                   fwhmValues = MALDIquant::metaData(isampled)$fwhm)
                  })

                  #// put everything together
                  fwhmdf       <- do.call("rbind", fwhmdf)
                  fwhmdf       <- fwhmdf[order(fwhmdf$peaks), ]


            }


      } else{

            if(!MALDIquant::isMassSpectrum(s)) {
                  stop("Input s is not a MassSpectrum object. See ?MALDIquant::MassSpectrum.\n")
            }

            #// detect peaks
            p           <- MALDIquant::detectPeaks(object = s, method = "SuperSmoother", SNR = 3)

            #// complain if peak is suspecious
            if(length(p@mass) == 0){
                  stop("The provided spectrum s does not contain any peaks. Consider providing another. \n")
            }

            if(length(p@mass) < 10){
                  warning("The provided spectrum s contains less than 10 peaks. Consider providing another. \n")
            }

            #// compute fwhm at every peak
            fwhmdf      <- .getFwhm(s, .samplep(p, peakSampling))

      }


      #// use super smoother to get a smooth curve
      sm           <- supsmu(x = fwhmdf$peaks, y = fwhmdf$fwhmValues, bass = 9)


      #// create the linear interpolation function
      fwhmFun      <- approxfun(x = sm$x, y = sm$y, rule = 2)

      if(plot)
      {
            r <- range(fwhmdf$peaks)
            qp <- seq(r[1], r[2])
            plot(x = fwhmdf$peaks, y = fwhmdf$fwhmValues,
                 main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
            lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

      }

      if(!is.null(savePlot))
      {
            # plot the the detecetd martix-pixels ontop of the optical image
            Cairo::Cairo(file = savePlot,
                         width = 3000, height = 2000, type = "svg",# units = "px",
                         dpi = 100)

            r <- range(fwhmFun$p)
            qp <- seq(r[1], r[2])
            plot(x = fwhmFun$p, y = fwhmFun$fwhmValues,
                 main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
            lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

            dev.off()
      }


      return(fwhm(fwhmInterpolator = fwhmFun,
                  peaks = fwhmdf$peaks,
                  fwhmVals = fwhmdf$fwhmValues))

}

#// internal function to sample detected peaks - for speed
# p is a MassPeaks object
.samplep      <- function(p, sampling) {

      if(sampling > length(p)){
            return(p)
      }

      idx         <- sample(seq(1, length(p)), sampling)
      orderedIdx  <- order(p@mass[idx])
      md          <- p@metaData # record metaData

      if(!is.null(p@metaData$fwhm)){
            p@metaData$fwhm <- p@metaData$fwhm[idx][orderedIdx]
      }

      MALDIquant::createMassPeaks(mass = p@mass[idx][orderedIdx],
                                  intensity = p@intensity[idx][orderedIdx],
                                  snr = p@snr[idx][orderedIdx],
                                  metaData = md)

}


#' Full width half Maximum of peaks
#'
#' A method to compute the fwhm for the peaks of a given spectrum.
#'
#' @param spectrum:    a spectrum of \code{MALDIquant::MassSpectrum} type.
#' @param peaks:  a centroided spectrum of \code{MALDIquant::MassPeaks} type.
#' @return Returns a named vector of fwhm values with the same length as the number of peaks within the spectrum.
#'
#'
#' @keywords internal
#'
#' @references Adopted from \url{https://gist.github.com/sgibb/3914291}.
#'
#'

.getFwhm              <- function(spectrum, peaks)
{

      #// complain if no peaks are found
      if(length(MALDIquant::mass(peaks)) < 1) {stop("peaks object is empty.\n")}

      #// if zero intensities encountered, add a small baseline value
      #zeros         <- MALDIquant::intensity(spectrum) == 0
      #if(any(zeros)) {MALDIquant::intensity(spectrum)[zeros] <- mean(MALDIquant::intensity(spectrum))}


      fwhmValues        <- sapply(MALDIquant::match.closest(x = MALDIquant::mass(peaks), MALDIquant::mass(spectrum)),
                                  .fwhm,
                                  spectrum = spectrum)

      fwhmdf            <- data.frame(peaks = MALDIquant::mass(peaks),
                                      fwhmValues = fwhmValues)


      return(fwhmdf)
}


## work horse for .getFwhm
.fwhm         <- function(spectrum, i)
{


      n             <- length(spectrum)
      left          <- ifelse(i <= 1, 1, i)
      right         <- ifelse(i >= n, n, i)


      intensities   <- MALDIquant::intensity(spectrum)
      peaksMasses   <- MALDIquant::mass(spectrum)



      hm            <- intensities[i]/2

      while (left > 1 && intensities[left] > hm)
      {
            left   <- left-1
      }

      while (right < n && intensities[right] > hm)
      {

            right  <- right+1
      }

      #// this to fix for the occasional constant-line artifacts at the beginning and end of spectra.
      # if((intensities[right - 1] == intensities[right]) | (intensities[left] == intensities[left + 1]))
      # {return(NA_real_)}

      ## interpolate x values
      xleft         <- approx(x=intensities[left:(left+1)],
                             y=peaksMasses[left:(left+1)], xout=hm)$y

      xright        <- approx(x=intensities[(right-1):right],
                             y=peaksMasses[(right-1):right], xout=hm)$y

      return(abs(xleft-xright))

}

