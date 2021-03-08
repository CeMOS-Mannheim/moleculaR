#' Approximation of fwhm(m/z)
#'
#' This function Creates an overall approximation of fwhm as a function of m/z for the entire dataset. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param sampling: Sample size (number of spectra), to reduce redundant overhead. 
#' @param plot: whether to plot. 
#' @param savePlot: whether to save a plot as svg. 
#' @param path: path to the saved plot "estimated fwhm(m/z)"
#' @return
#' Returns a linear interplotor function fwhm(m/z). 
#'
#' @export
#'
estimateFwhm  = function(x, sampling = 1000, plot = FALSE, savePlot = FALSE, path = file.path(getwd(), "estimated fwhm.svg")) {
       
       #// sample specta 
        sampling        = ifelse(length(x) < sampling, length(x), sampling)
       idx           = sample(seq_along(x), sampling, replace = FALSE)
       
       #// collect fwhm and corresponding masses. 
       collect       = lapply(x[idx], FUN = function(ispect) {
              
              data.frame(m = MALDIquant::mass(ispect),
                         fwhm = MALDIquant::metaData(ispect)$fwhm)
       })
       
       #// put everything together
       collect       = do.call("rbind", collect)
       collect       = collect[order(collect$m), ]
       
       
       #// use super smoother to get a smooth curve
       sm            = supsmu(x = collect$m, 
                              y = collect$fwhm)
       
       
       #// create the linear interpolation function
       fwhmFun       = approxfun(x = sm$x, y = sm$y, rule = 2)
       
       if(plot)
       {
              idx = sample(seq(1, nrow(collect)), size = 10000)
              r = range(collect$m)
              qp = seq(r[1], r[2])
              plot(x = collect$m[idx], y = collect$fwhm[idx], 
                   main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
              lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)
              
       }
       
       if(savePlot)
       {
              # plot the the detecetd martix-pixels ontop of the optical image
              Cairo::Cairo(file = path,
                           width = 3000, height = 2000, type = "svg",# units = "px",
                           dpi = 100) #"auto"
              idx = sample(seq(1, nrow(collect)), 10000)
              r = range(collect$m)
              qp = seq(r[1], r[2])
              plot(x = collect$m[idx], y = collect$fwhm[idx], 
                   main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm", 
                   cex = 2, cex.lab = 2, cex.axis = 2, cex.main = 3)
              lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 3)
              
              
              
              dev.off()
       }
       
       return(fwhmFun)
       
}



# ##// estimate fwhm across the mass axis
# plot(x = MALDIquant::mass(e$msDataPeaksNeg[[100]]), y = MALDIquant::metaData(e$msDataPeaksNeg[[100]])$fwhm, type = "l")
# tmp                  = supsmu(x = MALDIquant::mass(e$msDataPeaksNeg[[100]]), y = MALDIquant::metaData(e$msDataPeaksNeg[[100]])$fwhm)
# lines(x = tmp$x, y = tmp$y, col = "red")
# af = approxfun(x = tmp$x, y = tmp$y)
# lines(x = seq(100, 1200), y = test(seq(100, 1200)), col = "green")
# 
# test          = smooth.spline(x = MALDIquant::mass(e$msDataPeaksNeg[[100]]), y = MALDIquant::metaData(e$msDataPeaksNeg[[100]])$fwhm, df = 12)
# lines(x = e$msDataPeaksNeg[[100]]@mass, y = predict(test, x = MALDIquant::mass(e$msDataPeaksNeg[[100]]))$y, col = "green")
# 
# fwhm                 = sapply(e$msDataPeaksNeg, function(ispect) {
#        
#        ispect@metaData$fwhm
# })