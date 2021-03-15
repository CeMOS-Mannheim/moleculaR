#' Estimation of fwhm(m/z)
#'
#' This creates a linear interpolator function that approximats fwhm as a function
#' of m/z for the given single s.
#'
#' @param s: 	      a spectrum of \code{MALDIquant::Masss} type.
#' @param sampling:  an integer specifying how many detected p to consider for fwhm estimation.
#' @param plot:      whether to plot the result
#' @param savePlot:  either \code{NULL} or file path to save a plot as svg.
#' @param returnValues: if set to \code{TRUE}, returns additionally the calculated fwhm values.
#' @return
#' Returns a linear interplotor function fwhm(m/z). If \code{returnValues = TRUE} a list is returned
#' storing the fwhm(m/z) interpolating function in addition to the calculated fwhm values.
#'
#' @export
#'
#'
estimateFwhm      = function(s, sampling = 1000L, plot = FALSE, savePlot = NULL, returnValues = FALSE) {

      #// detect peaks
      p           = MALDIquant::detectPeaks(object = s, method = "SuperSmoother", SNR = 3)


      #// compute fwhm at every peak
      fwhmdf      = .getFwhm(s, .samplep(p, sampling))

      #// use super smoother to get a smooth curve
      sm           = supsmu(x = fwhmdf$p, y = fwhmdf$fwhmValues)


      #// create the linear interpolation function
      fwhmFun      = approxfun(x = sm$x, y = sm$y, rule = 2)

      if(plot)
     {
            r = range(fwhmdf$p)
            qp = seq(r[1], r[2])
            plot(x = fwhmdf$p, y = fwhmdf$fwhmValues,
                 main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
            lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

     }

       if(!is.null(savePlot))
       {
          # plot the the detecetd martix-pixels ontop of the optical image
          Cairo::Cairo(file = savePlot,
                       width = 3000, height = 2000, type = "svg",# units = "px",
                       dpi = 100)

           r = range(fwhmFun$p)
           qp = seq(r[1], r[2])
           plot(x = fwhmFun$p, y = fwhmFun$fwhmValues,
                main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
           lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

           dev.off()
       }

      if(returnValues) {
         return(list(fwhmValues = fwhmdf, fwhmFun = fwhmFun))
      } else{
         return(fwhmFun)
      }

}

#// internal function to sample detected peaks - for speed
.samplep      = function(p, sampling) {

   if(sampling > length(p)){
      return(p)
   }

   idx         = sample(seq(1, length(p)), sampling)
   orderedIdx  = order(p@mass[idx])

   MALDIquant::createMassPeaks(mass = p@mass[idx][orderedIdx],
                               intensity = p@intensity[idx][orderedIdx],
                               snr = p@snr[idx][orderedIdx])

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

.getFwhm              = function(spectrum, peaks)
{

   #// complain if no peaks are found
   if(length(MALDIquant::mass(peaks)) < 1) {stop("peaks object is empty.\n")}

   #// if zero intensities encountered, add a small baseline value
   #zeros         = MALDIquant::intensity(spectrum) == 0
   #if(any(zeros)) {MALDIquant::intensity(spectrum)[zeros] = mean(MALDIquant::intensity(spectrum))}


   fwhmValues                  = sapply(MALDIquant::match.closest(x = MALDIquant::mass(peaks), MALDIquant::mass(spectrum)),
                                        .fwhm,
                                        spectrum = spectrum)

   fwhmdf            = data.frame(peaks = MALDIquant::mass(peaks),
                                  fwhmValues = fwhmValues)


   return(fwhmdf)
}


## work horse for .getFwhm
.fwhm         = function(spectrum, i)
{


   n             = length(spectrum)
   left          = ifelse(i <= 1, 1, i)
   right         = ifelse(i >= n, n, i)


   intensities   = MALDIquant::intensity(spectrum)
   peaksMasses   = MALDIquant::mass(spectrum)



   hm            = intensities[i]/2

   while (left > 1 && intensities[left] > hm)
   {
      left   = left-1
   }

   while (right < n && intensities[right] > hm)
   {

      right  = right+1
   }

   #// this to fix for the occasional constant-line artifacts at the beginning and end of spectra.
   # if((intensities[right - 1] == intensities[right]) | (intensities[left] == intensities[left + 1]))
   # {return(NA_real_)}

   ## interpolate x values
   xleft         = approx(x=intensities[left:(left+1)],
                          y=peaksMasses[left:(left+1)], xout=hm)$y

   xright        = approx(x=intensities[(right-1):right],
                          y=peaksMasses[(right-1):right], xout=hm)$y

   return(abs(xleft-xright))

}


