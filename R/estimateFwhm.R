#' Estimation of fwhm as a function of m/z axis
#'
#' Estimate FWHM as a function of m/z from spectra, or use an approximate
#' instrument model when spectra are unavailable.
#'
#' @param s A numeric profile returned by [readSingleSpect()], a list/data frame
#'   with mz/intensity vectors, a list of those profiles, a native Cardinal
#'   object with known profile status, or an mz/fwhm table of precomputed widths
#'   (optional snr column). A numeric m/z vector can be paired with intensity.
#'   Existing MALDIquant spectra/lists remain temporarily supported. Omit or
#'   supply NULL to use [useDefaultFwhm()].
#' @param spectraSampling  an integer specifying how many spectra should be used for fwhm estimation
#' (the lower the faster) for a list or Cardinal experiment. Defaults to 10.
#' @param peakSampling  an integer specifying how many detected peaks (in each spectrum) to consider for fwhm estimation,
#' the lower the faster. Defaults to 1000.
#' @param numCores Positive worker count for native Cardinal processing. The
#'   default is serial; values above one use Snow workers. Legacy inputs retain
#'   their existing platform-specific behavior.
#' @param smoothingMethod a character, one of `c('loess','superSmoother')`, the method that is used for fitting a smoothing curve to
#' represent fwhm as a function of m/z. Defaults to `loess` which produced better fits but takes
#' longer to compute.
#' @param SNR an optional numeric, the minimum SNR level of the peaks that are used for fwhm estimation.
#' @param plot      whether to plot the result.
#' @param savePlot  either \code{NULL} or file path to save a plot as svg.
#' @param modality Instrument type used only when `s` is `NULL`: `"FTICR"`
#' (default), `"TOF-reflector"`, or `"TIMS-TOF"`. Aliases `"TOF"` and
#' `"timsTOF"` are accepted.
#' @param userMassRes Optional resolving power and reference m/z, supplied as
#' `c(massRes = 150000, mass = 400)`. Used only when `s` is `NULL`.
#' @param constant Logical; use a constant width instead of the instrument
#' model when `s` is `NULL`. Passed to [useDefaultFwhm()].
#' @param mz Optional numeric mass vector paired with intensity, instead of s.
#' @param intensity Numeric profile intensities paired with numeric s or mz.
#' @param peakIndices Optional unique sample indices for one numeric profile.
#'   Supplying these isolates peak-width calculations from automatic picking.
#' @param peakMethod Native Cardinal peakPick method; default diff. Other
#'   supported methods are sd, mad, quantile, filter and cwt.
#' @param BPPARAM Optional native BiocParallel backend; keep numCores at one when
#'   supplying it. Caller-started backends retain their lifecycle.
#'
#' @details
#' When `s` is omitted or `NULL`, [useDefaultFwhm()] supplies an approximate
#' model and a message reports the instrument and resolving power assumed.
#' Defaults are 123000 (FTICR), 10000 (TOF-reflector), and 40000 (TIMS-TOF),
#' all at m/z 249. FTICR scales the bundled nonlinear curve; TOF models use
#' constant resolving power. Measured widths are preferable when available.
#' Invalid spectra or failed estimation raise errors rather than using a default.
#' Sampling, smoothing and SNR settings apply only to supplied spectra.
#'
#' Native automatic picking uses Cardinal, which can select different peaks
#' from the historical MALDIquant SuperSmoother detector. The half-height linear
#' interpolation and the loess/SuperSmoother curve-fitting choices are preserved.
#' Sampling uses R's random stream; set.seed() makes selected spectra and peaks
#' reproducible. Only selected Cardinal profiles are retrieved. Profile masses
#' must be strictly increasing; select one run and complete queued processing
#' explicitly before fitting a Cardinal experiment.
#' Known centroided or unknown-status Cardinal input cannot supply measured
#' widths directly. Use a native mz/fwhm table when widths were measured elsewhere.
#' Invalid/unbounded calculated widths are discarded with a warning and recorded;
#' at least five distinct valid width samples are required for native fitting.
#'
#' If `s` is a list of centroided `MassPeaks` objects, the function assumes that FWHM values have been
#' (externally) pre-computed and attached to each centroided spectrum (`MassPeaks` object) via the
#' `metaData` slot (see `?MALDIquant::MassPeaks`) such that for every centroided spectrum `s[[i]]` the
#' FWHM values corresponding to all detected peaks `s[[i]]@mass` are stored in `s[[i]]@metaData$fwhm`.
#'
#' @return
#' An S3 `fwhm` object containing `fwhmInterpolator`, `peaks`, and `fwhmVals`.
#' For default models, widths are model predictions rather than measurements.
#' Native measured models also contain estimation metadata recording sampled
#' spectra, peak selection, detected/discarded counts and processing settings.
#'
#' @examples
#' fticr <- estimateFwhm()
#' timstof <- estimateFwhm(modality = "TIMS-TOF")
#' custom <- estimateFwhm(userMassRes = c(massRes = 150000, mass = 400))
#' getFwhm(custom, 400)
#' # With measured data: estimateFwhm(s = spectra)
#'
#' @export
#' @md
#'
estimateFwhm            <- function(s = NULL, spectraSampling = 10, peakSampling = 1000,
                                    numCores = 1, smoothingMethod = "loess",
                                    SNR = 5, plot = FALSE, savePlot = NULL,
                                    modality = "FTICR", userMassRes = NULL,
                                    constant = FALSE, mz = NULL, intensity = NULL,
                                    peakIndices = NULL, peakMethod = "diff", BPPARAM = NULL) {

      native <- .msiIsObject(s) || is.numeric(s) || !is.null(mz) || !is.null(intensity) ||
                .msiIsNumericSpectrum(s) ||
                (is.list(s) && length(s) > 0L && .msiIsNumericSpectrum(s[[1]]))
      if (native) {
            return(.estimateNativeFwhm(s, spectraSampling, peakSampling, numCores,
                   smoothingMethod, SNR, plot, savePlot, mz, intensity,
                   peakIndices, peakMethod, BPPARAM))
      }
      if (!is.null(peakIndices) || !is.null(BPPARAM)) {
            stop("peakIndices and BPPARAM require native measured spectra.")
      }

      if(is.null(s)){
            model <- useDefaultFwhm(userMassRes = userMassRes,
                                    modality = modality, constant = constant)
            instrument <- switch(modality, "TOF" = "TOF-reflector",
                                  "timsTOF" = "TIMS-TOF", modality)
            instrument <- match.arg(instrument, c("FTICR", "TOF-reflector", "TIMS-TOF"))
            referenceMass <- if(is.null(userMassRes)) 249 else {
                  if(is.null(names(userMassRes))) userMassRes[2] else userMassRes[["mass"]]
            }
            resolution <- referenceMass / getFwhm(model, referenceMass)
            message("No spectra supplied; using approximate ", instrument,
                    if(constant) " constant-width" else "",
                    " FWHM model (resolving power ", signif(resolution, 7),
                    " at m/z ", referenceMass, "). Measured widths are preferable when available.")
            .plotFwhmResult(model, plot, savePlot)
            return(model)
      }

      if(is.list(s)){

            if(MALDIquant::isMassSpectrumList(s)) {


                  #// sample spectra - here s is a list of MassPSpectrum objects
                  spectraSampling <- ifelse(length(s) < spectraSampling, length(s), spectraSampling)
                  idx         <- sample(x = seq_along(s), size = spectraSampling)
                  s           <- s[idx]

                  #// detect peaks - here p is a list of MassPeaks objects
                  p           <- MALDIquant::detectPeaks(object = s, method = "SuperSmoother", SNR = SNR)

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
                             "attached to s. Have you already applied peak binning? ",
                             "if so, please run fwhm estimation on raw data .. \n")
                  }

                  spectraSampling <- ifelse(length(s) < spectraSampling, length(s), spectraSampling)
                  idx         <- sample(x = seq_along(s), size = spectraSampling)
                  s           <- s[idx]

                  #// collect fwhm and corresponding masses.
                  fwhmdf       <- lapply(s, FUN = function(i) {

                        i <- .snrFilterp(i, SNR)

                        i <- .samplep(i, peakSampling)

                        data.frame(peaks = MALDIquant::mass(i),
                                   fwhmValues = MALDIquant::metaData(i)$fwhm)
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
            p           <- MALDIquant::detectPeaks(object = s, method = "SuperSmoother", SNR = SNR)

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

      #// create the linear interpolation function
      fwhmFun <- switch(smoothingMethod,
                        "superSmoother" = .superSmoother(x = fwhmdf$peaks, y = fwhmdf$fwhmValues),
                        "loess" = .loess(data = fwhmdf),
                        stop("smoothingMethod not recognized, accepted values = c('loess','superSmoother').\n"))

      model <- fwhm(fwhmInterpolator = fwhmFun,
                     peaks = fwhmdf$peaks, fwhmVals = fwhmdf$fwhmValues)
      .plotFwhmResult(model, plot, savePlot)
      return(model)

}

# Shared by measured estimates and default models.
.plotFwhmResult <- function(model, plot, savePlot){
      if(plot) graphics::plot(model)
      if(!is.null(savePlot)){
            grDevices::svg(filename = savePlot, width = 12, height = 8)
            on.exit(grDevices::dev.off(), add = TRUE)
            graphics::plot(model)
      }
      invisible(NULL)
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
            md$fwhm <- md$fwhm[idx][orderedIdx]
      }

      MALDIquant::createMassPeaks(mass = p@mass[idx][orderedIdx],
                                  intensity = p@intensity[idx][orderedIdx],
                                  snr = p@snr[idx][orderedIdx],
                                  metaData = md)

}

#// internal function to filter detected peaks according to their SNR
# p is a MassPeaks object
.snrFilterp      <- function(p, SNR) {

      if(all(is.na(p@snr))){
            warning("There is no SNR information attached to 's'. ",
                    "SNR peak filtration was skipped. \n")
            return(p)
      }

      if(min(p@snr, na.rm = TRUE) >= SNR){
            return(p)
      }

      idx         <- which(p@snr >= SNR)

      if(length(idx) == 0){
            warning(paste0("There are peaks with SNR >= ", SNR,
                           ". All available peaks will be used. \n"))
            return(p)
      }

      md          <- p@metaData # record metaData

      if(!is.null(md$fwhm)){
            md$fwhm <- md$fwhm[idx]
      }

      MALDIquant::createMassPeaks(mass = p@mass[idx],
                                  intensity = p@intensity[idx],
                                  snr = p@snr[idx],
                                  metaData = md)

}

#// compute the fwhm interpolating function by super smoother method
.superSmoother <- function(x, y){

      #// use super smoother to get a smooth curve
      sm           <- supsmu(x = x, y = y, bass = 10)

      #// create the linear interpolation function
      return(approxfun(x = sm$x, y = sm$y, rule = 2))
}

#// compute the fwhm interpolating function by loess method
.loess <- function(data){

      #// use loess to get a smooth curve
      lo           <- loess(fwhmValues ~ peaks, data,
                            control=loess.control(surface="direct"))

      #// create the linear interpolation function
      return(function(x) {predict(object = lo, newdata = x)})
}


#' Full width half Maximum of peaks
#'
#' A method to compute the fwhm for the peaks of a given spectrum.
#'
#' @param spectrum    a spectrum of \code{MALDIquant::MassSpectrum} type.
#' @param peaks  a centroided spectrum of \code{MALDIquant::MassPeaks} type.
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

