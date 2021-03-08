#' Gaussian weight as a function of distance and fwhm
#'
#' This function is used to estimate the appropriate Gaussian weighting to scale the intensities to be used with the lipid search method. 
#' The Gausiian is constructed from knowing the fwhm. 
#'
#' @param x: 	mass shift; the distance (in Da) from the theoretical mass that is being searched to the measured masses in a spectrum.
#' @param m: a numeric, the theoretical mass.  
#' @param fwhm: a numeric, the estimated fwhm as a function of m/z, see `?estimateFwhm`. 
#' @param plot: whether to plot. 
#' @param savePlot: whether to save a plot as svg. 
#' @param path: path to the saved plot `paste0(lipidSearch-m.svg)`
#' @return
#' Returns a numeric vector of the same length as x, representing the corresponding weights. 
#'
#' @export
#'
gaussWeight   = function(x, m, fwhm, ionMode = NA, plot = FALSE, savePlot = FALSE, path = file.path(getwd(), paste0("lipidSearch-", m,".svg"))) {
       
       #// create a Gaussian function that represents the theoretical peak based on the caluculated 
       #   fwhm function of the current dataset
      
       # fwhm = 2.355 * sigma
       s      = fwhm / 2.355
       
       w      = dnorm(x, mean = m, sd = s) / dnorm(m, mean = m, sd = s) # normalizatio is done to set the peak to 1 
       
       
       
       if(plot)
       {
              par(las = 2, mar = c(10, 5, 5, 5), cex = 3)
              gx     = seq(m-(10*s), m+(10*s), s/10)
              gy     = dnorm(gx, mean = m, sd = s) / dnorm(m, mean = m, sd = s)
               
              plot(y = gy, x = gx, type = "l", main = "Estimated Gaussian",  ylab = "Intensity (a.u.)", xlab = "",axes = FALSE, cex.lab = 1.5)
              
              ticks  = x
              
              if(is.na(ionMode)) {
                     axis(side = 1, at = x, col.ticks = "red", lwd = 2)
                     axis(side = 1, at = m, col.ticks = "black", lwd = 2)
              } else {
                     ticksNeg   = which(ionMode < 0)
                     ticksPos   = which(ionMode > 0)
                     
                     axis(side = 1, at = x[ticksNeg], col.ticks = "blue", lwd = 2)
                     axis(side = 1, at = x[ticksPos], col.ticks = "red", lwd = 2)
                     axis(side = 1, at = m, col.ticks = "black", lwd = 2)
              }
              
              title(xlab = "m/z (Da)", cex.lab = 1.5, line = 7)
              axis(side = 2)
              abline(v = m, col = "black", lty = "dashed")
              
              
              
              
       }
       
       if(savePlot)
       {
              # plot the the detecetd martix-pixels ontop of the optical image
              Cairo::Cairo(file = path,
                           width = 3000, height = 2000, type = "svg",# units = "px",
                           dpi = 100) #"auto"
              
              par(las = 2, mar = c(10, 5, 5, 5), cex = 3)
              gx     = seq(m-(10*s), m+(10*s), s/10)
              gy     = dnorm(gx, mean = m, sd = s) / dnorm(m, mean = m, sd = s)
              
              plot(y = gy, x = gx, type = "l", main = "Estimated Gaussian",  ylab = "Intensity (a.u.)", xlab = "",
                   axes = FALSE, cex.lab = 3, bty = "n")
              
              ticks  = x
              
              if(is.na(ionMode)) {
                     axis(side = 1, at = x, col.ticks = "red", lwd = 2, cex.lab = 3)
                     axis(side = 1, at = m, col.ticks = "black", lwd = 2, cex.lab = 3)
              } else {
                     ticksNeg   = which(ionMode < 0)
                     ticksPos   = which(ionMode > 0)
                     
                     axis(side = 1, at = x[ticksNeg], col.ticks = "blue", lwd = 2, cex.lab = 3)
                     axis(side = 1, at = x[ticksPos], col.ticks = "red", lwd = 2, cex.lab = 3)
                     axis(side = 1, at = m, col.ticks = "black", lwd = 2, cex.lab = 3)
              }
              
              title(xlab = "m/z (Da)", cex.lab = 1.5, line = 7)
              axis(side = 2, cex = 3)
              abline(v = m, col = "black", lty = "dashed")
              
              
              
              dev.off()
       }
       
       return(w)
       
}



