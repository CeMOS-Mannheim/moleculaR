#' Full width half Maximum of peaks
#'
#' A method to compute the fwhm for the peaks of a given spectrum.
#'
#' @param spectrum:    a spectrum of \code{MALDIquant::MassSpectrum} type.
#' @param peaks:  a centroided spectrum of \code{MALDIquant::MassPeaks} type.
#' @return Returns a named vector of fwhm values with the same length as the number of peaks within the spectrum.
#'
#' @export
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


## work horse
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
