#' Estimation of fwhm(m/z)
#'
#' This creates a linear interpolator function that approximats fwhm as a function
#' of m/z for the given single spectrum.
#'
#' @param spectrum: 	a spectrum of \code{MALDIquant::MassSpectrum} type.
#' @param peaks:        a centroided spectrum of \code{MALDIquant::MassPeaks} type.
#' @param plot:         whether to plot the result
#' @param savePlot:     either \code{NULL} or file path to save a plot as svg.
#' @return
#' Returns a linear interplotor function fwhm(m/z).
#'
#' @export
#'
#'
estimateFwhm        = function(spectrum, plot = FALSE, savePlot = NULL) {

      #// detect peaks
      peaks       = MALDIquant::detectPeaks(object = spectrum,
                                              method = "SuperSmoother",
                                              SNR = 3)
      #// compute fwhm at every peak
      fwhmdf      = .getFwhm(spectrum, peaks)

      #// use super smoother to get a smooth curve
      sm            = supsmu(x = fwhmdf$peaks,
                              y = fwhmdf$fwhmValues)


      #// create the linear interpolation function
      fwhmFun       = approxfun(x = sm$x, y = sm$y, rule = 2)

      if(plot)
     {
            r = range(fwhmdf$peaks)
            qp = seq(r[1], r[2])
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

           r = range(fwhmFun$peaks)
           qp = seq(r[1], r[2])
           plot(x = fwhmFun$peaks, y = fwhmFun$fwhmValues,
                main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
           lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

           dev.off()
       }

       return(fwhmFun)

}



