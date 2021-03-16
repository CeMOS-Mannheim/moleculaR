#' Calculate analyte probablity map 
#'
#' This function creates an analyte probability map from a given spatial point pattern. 
#'
#' @param x: 	       x values. 
#' @param y:         y values.
#' @param df:        degrees of freedom for the smoothing spline. 
#' @param adjRange:  whether to adjust the range of the second derivative to match
#' the range of y so it could be displayed along the original data. 
#' @param xQuery:    x values, at which to evaluate the fitted spline. 
#' A list .. 
#'
#' @export
#'
elbowPointSimple     = function(x, y, df = 7, adjRange = TRUE, 
                                xQuery = seq(range(x)[1], range(x)[2], 0.1)) {
       
       
       
       # fit a smoothing spline and fine the 2nd derivative. 
       smoothspl      = smooth.spline(x = x, y = y, df = df)
       
      
       deriv2        = predict(smoothspl, x = xQuery, deriv = 2)$y
       spl           = predict(smoothspl, x = xQuery)$y
       
       
       if(adjRange) { # adjust the range of the 2nd derivative to be able to plot it alongside x
              
              linMap <- function(x, a, b) approxfun(range(x), c(a, b))(x)
              
              deriv2   = linMap(deriv2, range(y)[1],  range(y)[2])
              
       }
       
       return(list(elbowPoint = xQuery[which.min(deriv2)],
                   x = xQuery,
                   deriv2 = deriv2, 
                   spl = spl))                                          
       
       
       
       
}



# # fit a smoothing spline and fine the 2nd derivative. 
# smoothsp      = with(bwdf, smooth.spline(x = bw, y = area, df = 7))
# 
# bwdf$deriv2   = predict(smoothsp, x = bwdf$bw, deriv = 2)$y
# bwdf$spline   = predict(smoothsp, x = bwdf$bw)$y
# 
# 
# linMap <- function(x, a, b) approxfun(range(x), c(a, b))(x)
# 
# bwdf$deriv2   = linMap(bwdf$deriv2, range(bwdf$area)[1],  range(bwdf$area)[2])
