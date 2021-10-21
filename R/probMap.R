#' Calculate molecular probablity maps
#'
#' This function creates molecular probability maps from a given spatial point pattern representation.
#'
#' @param sppMoi 	   The spatial point pattern.
#' @param bw        The Gaussian bandwidth used for Kernel density estimation. If a single numeric value is supplied
#' then it will be used for the computation of KDE. If a vector is supplied (in increasing order), then these values
#' will be passed to internal methods `spAutoCor` or `iterative` for computing the bandwidth based on cross validation,
#' see `?moleculaR::.bw.spAutoCorr` and `?moleculaR::.bw.iterative`.
#' @param bwMethod The method used for computing the Gaussian bandwidth, c("spAutoCor","iterative", "scott").
#' @param bwPlot    whether to plot the Gaussian bw selection procedure, ignored when \code{bwMethod = "scott"}.
#' @param csrIntensities How the intensities for the csrMoi model are generated, c("resample", "Poisson", "Gaussian").
#' @param sqrtTansform whether to apply square root transformation to the intensities of `sppMoi`. Set to 'TRUE' when
#' 'sppMoi' cotains a set of different analytes as in the case of collective projection maps.
#' @param control   An sppMoi object that is designated as the "control" or "null" alternative of \code{sppMoi}. If supplied the
#' csrMoi model will be generated from the intensities of this object.
#' @param csrMoi       Pre-computed wighted csrMoi model corresponding to the given sppMoi. This
#' could be useful when generating hotspot-bw curves.
#' @param pvalThreshold The p-value threshold to be used for the hypothesis testing.
#' @param pvalCorrection The method used for p-values correction, see '?p.adjust' for more details.
#' @param verbose   whether to output progress.
#' @return
#' #' An S3 object of type molProbMap with the following entries \cr
#' \itemize{
#'   \item bw: calculated Gaussian bandwidth.
#'   \item sppMoi: the input spp with intensities sqrt-transformed if 'sqrtTransform == TRUE'.
#'   \item csrMoi: calculated complete spatial randomness model for the input sppMoi.
#'   \item rhoMoi: density image of the input point pattern sppMoi.
#'   \item rhoCsr: density image of the computed csrMoi pattern.
#'   \item hotspotpp: the remaining points (type 'ppp') which lie within the hotspot mask as a spatial point pattern.
#'   \item hotspotIm: the hotspot image of type 'im'.
#'   \item hotspotMask: the hotspot mask of type 'owin'.
#'   \item coldspotpp: the remaining points (type 'ppp') which lie within the coldspot mask as a spatial point pattern.
#'   \item coldspotIm: the coldspot image of type 'im'.
#'   \item coldspotMask: the hotspot mask of type 'owin'.
#' }
#'
#'
#' @export
#'
probMap                     <- function(sppMoi,  control = NULL,
                                       bw = seq(1, 10, 1), bwMethod = "spAutoCor",
                                       csrIntensities = "resample", sqrtTansform = FALSE,
                                       csrMoi = NULL, pvalThreshold = 0.05, pvalCorrection = "bonferroni",
                                       bwPlot = FALSE, verbose = TRUE) {



      if(!("analytePointPattern" %in% class(sppMoi))){
         stop("sppMoi must be of 'analytePointPattern' and 'ppp' class. \n")
      }

       if(is.null(csrMoi)) {


             if(sqrtTansform){
                   sppMoi$marks$intensity <- sqrt(sppMoi$marks$intensity)
             }

             if(is.null(control)){ # if control is not supplied generate the csrMoi model from sppMoi

                   csrMoi              <- spatstat::rpoint(n = sppMoi$n, win = sppMoi$window)


                   csrMoi$marks   <- switch(csrIntensities,
                                          "resample" = {
                                                data.frame(intensity = sample(sppMoi$marks$intensity))
                                          },
                                          "Poisson" = {
                                                data.frame(intensity = rpois(csrMoi$n, mean(sppMoi$marks$intensity)))
                                          },
                                          "Gaussian"= {
                                                data.frame(intensity = rnorm(csrMoi$n, mean(sppMoi$marks$intensity),
                                                                              sd(sppMoi$marks$intensity)))
                                          })


                   csrMoi <- analytePointPattern(csrMoi, sppMoi$mzVals)

             } else {

                   if(!(spatstat::is.ppp(control))) {stop("Supplied control is not an ppp object. \n")}


                   if(!("analytePointPattern" %in% class(control))){
                      stop("control must be of 'analytePointPattern' and 'ppp' class. \n")
                   }

                   if(sqrtTansform){
                         sppMoi$marks$intensity <- sqrt(sppMoi$marks$intensity)
                         control$marks$intensity <- sqrt(control$marks$intensity)
                   }

                   csrMoi              <- spatstat::rpoint(n = sppMoi$n, win = sppMoi$window)



                   if(is.null(control$marks$intensity)){
                         stop("The supplied control does not contain intensity info in its marks.\n")
                   }

                   csrMoi$marks   <- switch(csrIntensities,
                                        "resample" = {
                                              if(control$n >= csrMoi$n){
                                                    data.frame(intensity = sample(control$marks$intensity,
                                                                                  size = csrMoi$n,
                                                                                  replace = FALSE))
                                              } else {
                                                    data.frame(intensity = sample(control$marks$intensity,
                                                                                  size = csrMoi$n,
                                                                                  replace = TRUE))
                                              }

                                        },
                                        "Poisson" = {
                                              data.frame(intensity = rpois(csrMoi$n, mean(control$marks$intensity)))
                                        },
                                        "Gaussian"= {
                                              data.frame(intensity = rnorm(csrMoi$n, mean(control$marks$intensity),
                                                                           sd(control$marks$intensity)))
                                        })

                   csrMoi <- analytePointPattern(csrMoi, sppMoi$mzVals)



              }



       } else { # if csrMoi is supplied check if intensity-marks are supplied

          if(!("analytePointPattern" %in% class(csrMoi))){
             stop("csrMoi must be of 'analytePointPattern' and 'ppp' class. \n")
          }

          if(is.null(csrMoi$marks$intensity) )
                stop("The supplied csrMoi model does not contain intensity info in its marks.\n")

       }


       # compute kernel density bandwidth
       if(length(bw) > 1)
       {
             if(any(diff(bw) < 0))
                   stop("supplied bw must be either a single numeric of a vector of numerics of increasing order.\n")

              bw        <- switch(bwMethod,
                                 spAutoCor = {
                                         .bw.spAutoCorr(sppMoi = sppMoi,
                                                           plot = bwPlot,
                                                           bw = bw,
                                                           verbose = verbose)
                                 },
                                 iterative = {
                                         .bw.iterative(sppMoi = sppMoi,
                                                       bw = bw,
                                                       csrMoi = csrMoi, plot = bwPlot,
                                                       verbose = verbose)
                                 },
                                 scott = {
                                         .bw.scott.iso2(sppMoi)
                                 },
                                 stop("wrong bwMethd specified, must be one of c('spAutoCor','iterative','scott')"))



       }



      csrMoiw       <- csrMoi$marks$intensity
      sppMoiw       <- sppMoi$marks$intensity


       win <- spatstat::as.mask(sppMoi$window,
                              dimyx=c(diff(sppMoi$window$yrange),diff(sppMoi$window$xrange)))


       # create a density map for csrMoi
       rhoCsr           <- spatstat::density.ppp(x = csrMoi, sigma = bw,
                                                weights = csrMoiw, W = win, positive = TRUE)


       # scale such that sum{pixels} <= 1 i.e. a probability density function
       rhoCsr           <- rhoCsr/sum(rhoCsr)

       # create a density map for the image
       rhoMoi           <- spatstat::density.ppp(x = sppMoi, sigma = bw,
                                                 weights = sppMoiw, W = win, positive = TRUE)

       # scale such that sum{pixels} <= 1 i.e. a probability density function
       rhoMoi           <- rhoMoi/sum(rhoMoi)




       ## __ statistical testing and pvalue correction __ #

       # null hypothesis
       mucsrMoi            <- mean(rhoCsr, na.rm = TRUE)
       sigmacsrMoi         <- sd(rhoCsr, na.rm = TRUE)

       # convert to data.frame
       rhoMoidf <- spatstat::as.data.frame.im(x = rhoMoi)
       pvalsLwr <- rhoMoidf
       pvalsUpr <- rhoMoidf

       # generate p-values - lower tail & upper tail
       pvalsLwr$value <- pnorm(rhoMoidf$value, mean = mucsrMoi, sd = sigmacsrMoi, lower.tail = TRUE)
       pvalsUpr$value <- pnorm(rhoMoidf$value, mean = mucsrMoi, sd = sigmacsrMoi, lower.tail = FALSE)

       # correction
       pvalsLwr$value <- p.adjust(p = pvalsLwr$value, method = pvalCorrection)
       pvalsUpr$value <- p.adjust(p = pvalsUpr$value, method = pvalCorrection)

       # convert back to image
       pvalsLwr <- spatstat::as.im.data.frame(pvalsLwr)
       pvalsUpr <- spatstat::as.im.data.frame(pvalsUpr)


       ## __ hotspot __ ##

       hotspotIm        <- spatstat::eval.im(rhoMoi * (pvalsUpr <= pvalThreshold))

       # filter out points lying outside the computed hotspot
       tmpIm            <- hotspotIm
       tmpIm$v[which(tmpIm$v == 0, arr.ind = TRUE)] <- NA # manually set zero pixels to NA to remove them from mask
       hotspotMask      <- spatstat::as.owin(tmpIm)


       hotspotpp        <- spatstat::ppp(x = sppMoi$x,
                                        y = sppMoi$y,
                                        window = spatstat::as.polygonal(hotspotMask),
                                        marks = sppMoi$marks,
                                        checkdup = FALSE)

       hotspotpp        <- analytePointPattern(spp = hotspotpp, mzVals = sppMoi$mzVals) # only for conformity with S3 class


       ## __ coldspot __ ##

       coldspotIm        <- spatstat::eval.im(rhoMoi * (pvalsLwr <= pvalThreshold))

       # filter out points lying outside the computed coldspot
       tmpIm             <- coldspotIm
       tmpIm$v[which(tmpIm$v == 0, arr.ind = TRUE)] <- NA # manually set zero pixels to NA to remove them from mask
       coldspotMask      <- spatstat::as.owin(tmpIm)


       coldspotpp        <- spatstat::ppp(x = sppMoi$x,
                                         y = sppMoi$y,
                                         window = spatstat::as.polygonal(coldspotMask),
                                         marks = sppMoi$marks,
                                         checkdup = FALSE)

       coldspotpp        <- analytePointPattern(spp = coldspotpp, mzVals = sppMoi$mzVals) # only for conformity with S3 class


       return(molProbMap(bw = bw,		# calculated Gaussian bandwidth
                   sppMoi = sppMoi,       # the input spp with intensiteis sqrt-transfomed if sqrtTransform = TRUE
                   csrMoi = csrMoi,       # the created CSR model of the input spp
                   rhoMoi = rhoMoi,	# density image of the input point pattern
                   rhoCsr = rhoCsr,	# density image of the computed csrMoi pattern
                   hotspotpp = hotspotpp, # the remaining points which lie within the hotspot mask
                   hotspotIm = hotspotIm, # the hotspot image
                   hotspotMask = hotspotMask, # the hotspot mask of type owin.mask
                   coldspotpp = coldspotpp, # the remaining points which lie within the coldspot mask
                   coldspotIm = coldspotIm, # the coldspot image
                   coldspotMask = coldspotMask # the coldspot mask of type owin.mask
                   ))




}


#' Calculate Gaussian bandwidth based on spatial autocorrelation - `spAutoCor`
#'
#' This function computes the Gaussian bandwidth based on sensitivity analysis of spatial autocorrlation.
#' This is used internally in `moleculaR::probMap`.
#'
#' @param sppMoi: 	   The spatial point pattern, object of type `ppp`.
#' @param bw:        bandwidth steps at which to compute density and the corresponding Moran's I statistic.
#' @param plot:      whether to plot the result.
#' @param verbose:   whether to output progress.
#' @return
#' A numeric, the estimated Gaussian bandwidth.
#'
#' @export
#' @keywords internal

.bw.spAutoCorr <- function(sppMoi, bw = seq(1, 10, 1),
                           plot = FALSE, verbose = FALSE) {


      #// create a dataframe to hold the results
      bwdf        <- data.frame(bw = bw, moransI = numeric(length(bw)))
      sppMoiw        <- sppMoi$marks$intensity
      win         <- spatstat::as.mask(sppMoi$window, dimyx=c(diff(sppMoi$window$yrange),diff(sppMoi$window$xrange)))


      if(is.null(sppMoiw)){
            stop("sppMoi does not have intensity weights")
      }

      #// to show progress
      if(verbose)
            pb    <- utils::txtProgressBar(min = min(bw), max = max(bw), width = 20, style = 3)



      bwdf        <- lapply(bw, function(bwi){

            if(verbose)
                  utils::setTxtProgressBar(pb, bwi)




            # create a density map for the image
            rhoMoi      <- spatstat::density.ppp(x = sppMoi, sigma = bwi,
                                                 weights = sppMoiw,
                                                 W = win)

            moransppMoi    <- raster::Moran(raster::raster(rhoMoi), w = matrix(c(1,1,1,1,0,1,1,1,1), nrow=3))


            return(data.frame(bw = bwi, moransI = moransppMoi))

      })

      if(verbose)
            close(pb)

      bwdf <- do.call("rbind", bwdf)



      #// compute the elbow point
      ep                <- .inflictPoint(x = bwdf$bw, y = bwdf$moransI, plot = plot)


      if(verbose)
            cat("done.\n")

      return(c(bw = ep$inflictPoint))

}


#' Calculate Gaussian Bandwidth - `iterative`
#'
#' This function calculates the Gaussian bandwidth to be used for analyte probability maps by
#' finding the curve infliction point for a curve that represents hotspot area as a function
#' of bandwidth. This is used internally in `moleculaR::probMap`.
#'
#' @param sppMoi: 	The spatial point pattern.
#' @param bw:        A vector, The gaussian band width pool used for Kernel density estimation.
#' @param csrMoi:       Pre-computed wighted csrMoi model corresponding to the given sppMoi. This
#' speeds up the computation and is mostly used for troubleshooting.
#' @param pvalThreshold: The p-value threshold to be used for the hypothesis testing.
#' @param pvalCorrection: The method used for p-values correction, see '?p.adjust' for more details.
#' @param plot:      Whether to plot the hotspot area as function of bandwidth.
#' @param plotEach:  whether to plot the resulting significance area of each bandwidth iteration.
#' @param verbose:   Whether to show progress.
#' @return
#' A list ..
#'
#' @export
#' @keywords internal
#'

.bw.iterative                 <- function(sppMoi,  bw = seq(1, 10, 1), csrMoi = NULL,
                                         pvalThreshold = 0.05, pvalCorrection = "bonferroni",
                                         plot = FALSE, plotEach = FALSE, verbose = FALSE) {



       if(is.null(csrMoi)) {

             ## craete a complete spatial randomness point pattern with the same number of points and window ----
             csrMoi        <- spatstat::rpoint(n = sppMoi$n, win = sppMoi$window)

             csrMoi$marks  <- data.frame(intensity = sample(sppMoi$marks$intensity))

       } else { # if csrMoi is supplied check if intensity-marks are supplied

             if(is.null(csrMoi$marks$intensity))
                   stop("The supplied csrMoi model does not contain intensity info in its marks.\n")

       }



      csrMoiw       <- csrMoi$marks$intensity
      sppMoiw       <- sppMoi$marks$intensity


       #// create a dataframe to hold the results
       bwdf             <- data.frame(bw = bw, area = NA_real_)
       win              <- spatstat::as.mask(sppMoi$window, dimyx=c(diff(sppMoi$window$yrange),diff(sppMoi$window$xrange)))


       #// to show progress
       if(verbose)
             pb               <- utils::txtProgressBar(min = min(bw), max = max(bw), width = 20, style = 3)



       bwdf$area              <- sapply(bw, function(bwi){

             if(verbose)
                   utils::setTxtProgressBar(pb, bwi)



             # create a density map for csrMoi
             rhoCsr               <- spatstat::density.ppp(x = csrMoi, sigma = bwi,
                                                          weights = csrMoiw,
                                                          W = win)

             # scale such that sum{pixels} <= 1 i.e. a probability density function
             rhoCsr           <- rhoCsr/sum(rhoCsr)


             # create a density map for the image
             rhoMoi               <- spatstat::density.ppp(x = sppMoi, sigma = bwi,
                                                          weights = sppMoiw,
                                                          W = win)

             # scale such that sum{pixels} <= 1 i.e. a probability density function
             rhoMoi           <- rhoMoi/sum(rhoMoi)


             ## __ statistical testing and pvalue correction __ #

             # null hypothesis
             mucsrMoi            <- mean(rhoCsr, na.rm = TRUE)
             sigmacsrMoi         <- sd(rhoCsr, na.rm = TRUE)

             # convert to data.frame
             rhoMoidf <- spatstat::as.data.frame.im(x = rhoMoi)
             pvalsUpr <- rhoMoidf

             # generate p-values -  upper tail
             pvalsUpr$value <- pnorm(rhoMoidf$value, mean = mucsrMoi, sd = sigmacsrMoi, lower.tail = FALSE)

             # correction
             pvalsUpr$value <- p.adjust(p = pvalsUpr$value, method = pvalCorrection)

             # convert back to image
             pvalsUpr <- spatstat::as.im.data.frame(pvalsUpr)


             ## __ hotspot __ ##

             hotspotIm        <- spatstat::eval.im(rhoMoi * (pvalsUpr <= pvalThreshold))

             hotspotIm[hotspotIm == 0] <- NA # set zeros to NA to create a window
             hotspotWin       <- spatstat::as.polygonal(spatstat::as.owin(hotspotIm))

             if(plotEach){
                   par(mfrow = c(1,1))
                   spatstat::plot.owin(sppMoi$window, ylim = rev(sppMoi$window$yrange), add = FALSE, main = paste0("BW = ",bwi))
                   spatstat::plot.owin(hotspotWin, col = rgb(0,1,0,1), add = TRUE)
             }


             return(spatstat::area.owin(hotspotWin) / spatstat::area.owin(sppMoi$window))

       })

       if(verbose)
             close(pb)


       #// compute the elbow point
       if(all(bwdf$area == 0)) { # no hotspot detected

              cat("no hotspot detected for all bw. Setting arbitrary bw .. \n")
              # return(list(bwdf = bwdf,
              #             inflictPointData = list(inflictPoint = 3))) # arbitrary bw

             return(c(bw = 3)) # arbitrary bw
       }

       ep                          <- .inflictPoint(x = bwdf$bw, y = bwdf$area,
                                                   plot = plot, ylabel = "Relative Area")


       if(verbose)
             cat("done.\n")

       # return(list(bwdf = bwdf,
       #             inflictPointData = ep))

       return(c(bw = ep$inflictPoint))



}

#' Inflection point of a curve
#'
#' This function calculates the infliction point of a curve based on the 2nd derivative.
#' This is used internally in `moleculaR::.calcGaussBW`.
#'
#' @param x: 	   x values.
#' @param y:         y values.
#' @param df:        degrees of freedom for the smoothing spline.
#' @param adjRange:  whether to adjust the range of the second derivative to match
#' the range of y so it could be displayed along the original data.
#' @param plot:   whether to plot the result.
#' @param xQuery:    x values, at which to evaluate the fitted spline.
#' A list ..
#'
#' @export
#' @keywords internal
#'
.inflictPoint     <- function(x, y, df = 7, adjRange = TRUE, plot = FALSE,
                             xQuery = seq(range(x)[1], range(x)[2], 0.1),
                             ylabel = "Moran's I") {



      # fit a smoothing spline and fine the 2nd derivative.
      smoothspl      <- smooth.spline(x = x, y = y, df = df)


      #deriv1        <- predict(smoothspl, x = xQuery, deriv = 1)$y
      #deriv2        <- predict(smoothspl, x = xQuery, deriv = 2)$y
      spl           <- predict(smoothspl, x = xQuery)$y


      if(adjRange) { # adjust the range of the  derivative to be able to plot it alongside x

            linMap <- function(i, a, b) approxfun(range(i), c(a, b))(i)

            #deriv1   <- linMap(deriv1, range(spl)[1],  range(spl)[2])
            #deriv2   <- linMap(deriv2, range(spl)[1],  range(spl)[2])

      }

      mxCurve <- .getMaxCurve(xQuery, spl)
      #mxCurve <- .getMaxCurve(xQuery, deriv1)
      #mxCurve <- .getMaxCurve(xQuery, deriv2)



      if(plot){
            plot(x = xQuery, y = spl, type = "b", main = "Maximum Curvature",
                 xlab = "bw",
                 ylab = ylabel)


            #lines(x = xQuery, y = deriv1, lty = "solid" ,  col = rgb(0,0,1,0.5), lwd = 2)
            #lines(x = xQuery, y = deriv2, lty = "solid" ,  col = rgb(0,1,0,0.5), lwd = 2)

            abline(v = mxCurve, col = "chocolate1", lty ="dashed", lwd = 2)

            legend("right", bty = "n",
                   legend = c(paste0(ylabel," (bw)"),
                              #"1st-derivative",
                              #"2nd-derivative",
                              paste0("infliction=", round(mxCurve, 2))
                   ),
                   col = c("black",
                           #"blue",
                           #"green",
                           "chocolate1"
                   ),
                   lty = c("solid",
                           #"solid",
                           #"solid",
                           "dashed"
                   ))
      }

      return(list(inflictPoint = mxCurve,
                  x = xQuery,
                  #deriv1 = deriv1,
                  #deriv2 = deriv2,
                  spl = spl))




}



#'  Maximum curvature
#'
#' Finds the maximum curvature (elbow point) of vector of values representing a curve by finding the maximum distance from the
#' curve to the line drawn between the peak and tail of that curve.
#'
#' @param x:    a numeric vector.
#' @param plot: a logical to plot the result for validation purposes. Defaults to FALSE.
#'
#'
#' @return Returns the maximum curvature threshold (elbow point) of `x`.
#'
#' @export
#' @keywords internal
#'
#' @author Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
#'
#' @references \url{https://en.wikipedia.org/wiki/Unimodal_thresholding}
#' @references \url{https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line}
#'
#'

.getMaxCurve                      <- function(x, y, plot = FALSE)
{

      if(length(x) != length(y))
            stop("x and y have different lengths.")

      #// find the coordinates of the max bin (peak) and the last bin (tail)
      hpIdx                        <- which.max(y) # high point in y
      lpIdx                        <- which.min(y) # low point in y

      hp                          <- c(x = x[hpIdx], y = y[hpIdx])
      lp                          <- c(x = x[lpIdx], y = y[lpIdx])

      allPoints                   <- data.frame(x = x, y = y)

      searchSpace                 <- unname(unlist(apply(allPoints,
                                                        1,
                                                        FUN = function(x, p1 = hp, p2 = lp)
                                                        {
                                                              # ref: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line #
                                                              # suppose you have a line defined by two points P1 and P2, the the
                                                              # distance from point x to the line is defined (in 2D) by

                                                              abs(((p2[2] - p1[2]) * x[1]) - ((p2[1] - p1[1]) * x[2]) + (p2[1] * p1[2]) - (p2[2] * p1[1]) /
                                                                        sqrt(sum((p2 - p1) ** 2)))

                                                        })))

      if(lpIdx < hpIdx){

            if(lpIdx != 1)
                  searchSpace[1:lpIdx] <- NA
            if(hpIdx != length(x))
                  searchSpace[hpIdx:length(x)] <- NA

      } else {

            if(hpIdx != 1)
                  searchSpace[1:hpIdx] <- NA
            if(lpIdx != length(x))
                  searchSpace[lpIdx:length(x)] <- NA

      }

      maxCurveIdx                 <- which.max(searchSpace)

      r                            <- x[maxCurveIdx]


      if(plot)
      {
            plot(x = x, y = y, type = "l", main = "Maximum Curvature",
                 xlab = "x", ylab = "y")
            points(rbind(hp, lp, data.frame(x = r, y = y[maxCurveIdx])),
                   col = c("red", "red", "green"), cex = 1, pch = 19)
            lines(x = rbind(hp, lp), col = "red", lty = 2)
      }


      return(r)
}







#' Calculate Gaussian bandwideth - `scott`
#'
#' This function Calculate Gaussian bandwideth based on Scott's rule of thumb.
#' Note that this function is already available in recent \code{spatstat} versions (see \code{?spatstat::.bw.scott}).
#' This is used internally in `moleculaR::probMap`.
#'
#' @param X: 	   A point pattern (object of class \code{ppp}).
#' @param isotropic:       Logical value indicating whether to compute a single bandwidth
#' for an isotropic Gaussian kernel (isotropic=TRUE) or separate bandwidths for each
#' coordinate axis (isotropic=FALSE, the default).
#' @return
#' A numeric, estimated Gaussian bandwidth.

#// the following function is already available in recent spatstat versions.
.bw.scott2 <- function(X, isotropic=FALSE, d = NULL) {

       if(is.null(d)) { d <- spatstat::spatdim(X) } else spatstat.utils::check.1.integer(d)
       nX <- spatstat::npoints(X)
       cX <- spatstat::coords(X, spatial=TRUE, temporal=FALSE, local=FALSE)
       sdX <- apply(cX, 2, sd)
       if(isotropic) {
              #' geometric mean
              sdX <- exp(mean(log(pmax(sdX, .Machine$double.eps))))
       }
       b <- sdX * nX^(-1/(d+4))
       names(b) <- if(isotropic) "sigma" else paste0("sigma.", colnames(cX))
       return(b)
}

.bw.scott.iso2 <- function(X) { .bw.scott2(X, isotropic=TRUE) }



