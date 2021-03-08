#' Remove Peaks that are <= SNR
#'
#' This function removes peaks that have SNR <= `snr`. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param snr: The snr value. Defaults to 5. 
#' @return
#' A filtered list of `MassPeaks` objects. 
#'
#' @export
#'
peakSnrThresholding  = function(x, snr = 5) {
       
       
       
              x      = lapply(x, FUN = function(i) {


              idx           = which(i@snr >= 5)
              mt            = MALDIquant::metaData(i)
              mt$fwhm       = mt$fwhm[idx]

              MALDIquant::createMassPeaks(mass = i@mass[idx],
                                          intensity = i@intensity[idx],
                                          snr = i@snr[idx],
                                          metaData = mt)

       })
       
       
       
       
       return(x)
       
}

