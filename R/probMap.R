#' Calculate molecular probablity maps
#'
#' This function creates molecular probability maps from a given spatial point pattern representation.
#'
#' @param spp: 	   The spatial point pattern.
#' @param win:       The window object of type `spatstat::owin`.
#' @param weighted:  Whether the intensities are inlcude into the computation.
#' @param bw:        The Gaussian bandwidth used for Kernel density estimation.
#' @param bwMethod: The method used for computing the Gaussian bandwidth, c("spAutoCor","iterative", "scott").
#' @param dists:     for `bwMethod="spAutoCor"`, these are neighborhood distances at which to compute the Moran's I statisitc.
#' for `bwMethod="iterative"` these are the Gaussian bandwidth upon which to iterate.
#' @param numSim:    the number of Monte Carlo simulations to be done for Moran's I, see ?spdep::moran.mc.
#'   (only relevant with `bwMethod="spAutoCor"`).
#' @param approximate:    for `bwMethod="spAutoCor"`, when `TRUE`only points within a circle centered at the intensity-wighted
#' centroid (center of gravity) and radius `approxRadius` are considered. This is used to speed up the caluclation.
#' @param approxRadius: for `bwMethod="spAutoCor"`, the radius in units (usually pixels) for the approximation circlular region.By default,
#' this equates to the radius that produces a circular area which is half of the area of the spp window.
#' @param bwPlot:    whether to plot the Gaussian bw selection procedure, ignored when \code{bwMethod = "scott"}.
#' @param csrIntensities: How the intensities for the CSR model are generated, c("resample", "Poisson", "Gaussian").
#' @param control:   An spp object that is designated as the "control" or "null" alternative of \code{spp}. If supplied the
#' CSR model will be generated from the intensities of this object.
#' @param csr:       Pre-computed wighted csr model corresponding to the given spp. This
#' could be useful when generating generating hotspot-bw curves.
#' @param verbose:   whether to output progress.
#' @return
#' A list of
#' cutoff: the chosen threshold above which it is a hotspot,
#' GausBW: calculated Gaussian bandwidth,
#' csrspp: calculated complete spatial randomness model for the input spp,
#' denspp: density image of the input point pattern spp,
#' denCsr: density image of the computed csr pattern,
#' hotspotpp: the remaining points which lie within the hotspot mask as a spatial point pattern,
#' hotspotIm: the hotspot image,
#' hotspotMask: the hotspot mask of type owin.mask,
#' nonHotspotMask: the non-hotspot areas, to be used to mask the non hotspot areas.

#'
#' @export
#'
probMap                     = function(spp, win, weighted = TRUE, control = NULL,
                                       bw = NULL, bwMethod = "spAutoCor", dists = seq(1, 5, 0.5),
                                       numSim = 99, approximate = FALSE, approxRadius = NULL,
                                       csrIntensities = "resample", csr = NULL, bwPlot = FALSE,
                                       verbose = TRUE) {




       if(is.null(csr)) {

              ## craete a complete spatial randomness point pattern with the same number of points and window ----

             if(is.null(control)){ # if control is not supplied generate the csr model from spp

                   csr              = spatstat::rpoint(n = spp$n, win = win)

                   if(weighted){
                         csr$marks   = switch(csrIntensities,
                                                "resample" = {
                                                      data.frame(intensity = sample(spp$marks$intensity))
                                                },
                                                "Poisson" = {
                                                      data.frame(intensity = rpois(csr$n, mean(spp$marks$intensity)))
                                                },
                                                "Gaussian"= {
                                                      data.frame(intensity = rnorm(csr$n, mean(spp$marks$intensity),
                                                                                    sd(spp$marks$intensity)))
                                                })
                   }

             } else {

                   if(!(spatstat::is.ppp(control))) {stop("Supplied control is not an ppp object. \n")}
                   csr               = spatstat::rpoispp(lambda = spatstat::intensity(control), win = win)
                   #csr              = spatstat::rpoint(n = spp$n, win = win)


                   if(weighted){
                         if(is.null(control$marks$intensity) & weighted == TRUE){
                               stop("The supplied control does not contain intensity info in its marks.\n")
                               }
                         csr$marks   = switch(csrIntensities,
                                              "resample" = {
                                                    if(control$n >= csr$n){
                                                          data.frame(intensity = sample(control$marks$intensity,
                                                                                        size = csr$n,
                                                                                        replace = FALSE))
                                                    } else {
                                                          data.frame(intensity = sample(control$marks$intensity,
                                                                                        size = csr$n,
                                                                                        replace = TRUE))
                                                    }

                                              },
                                              "Poisson" = {
                                                    data.frame(intensity = rpois(csr$n, mean(spp$marks$intensity)))
                                              },
                                              "Gaussian"= {
                                                    data.frame(intensity = rnorm(csr$n, mean(spp$marks$intensity),
                                                                                 sd(spp$marks$intensity)))
                                              })

                   }
              }



       } else { # if csr is supplied check if intensity-marks are supplied if 'weighted == TRUE'

             if(is.null(csr$marks$intensity) & weighted == TRUE)
                   stop("The supplied csr model does not contain intensity info in its marks. Consider switching 'weighted' to FALSE.\n")

       }


       # compute kernel density bandwidth
       if(is.null(bw))
       {

              bw        = switch(bwMethod,
                                 spAutoCor = {
                                         .bw.spAutoCorr(spp = spp, win = win,
                                                           plot = bwPlot,
                                                           dists = dists,
                                                           numSim = numSim,
                                                           approximate = approximate,
                                                           approxRadius = approxRadius,
                                                           verbose = verbose)
                                 },
                                 iterative = {
                                         .bw.iterative(spp = spp, win = win, weighted = weighted,
                                                     csr = csr, plot = bwPlot,
                                                     verbose = verbose)
                                 },
                                 scott = {
                                         .bw.scott.iso2(spp)
                                 })



       }


      if(weighted){ # if the model is intensity-weighted assign the weights
            csrw       = csr$marks$intensity
            sppw       = spp$marks$intensity
      } else { # otherwise switch weights off
            csrw       = NULL
            sppw       = NULL
      }

       # create a density map for csr
       denCsr               = spatstat::density.ppp(x = csr, sigma = bw,
                                                    weights = csrw,
                                                    W = spatstat::as.mask(win,
                                                                          dimyx=c(diff(win$yrange) + 1,
                                                                                  diff(win$xrange) + 1)))
       # scale such that sum{pixels} <= 1 i.e. a probability density function
       denCsr           = denCsr/sum(denCsr)

       # create a density map for the image
       denspp               = spatstat::density.ppp(x = spp, sigma = bw,
                                                    weights = sppw,
                                                    W = spatstat::as.mask(win,
                                                                          dimyx=c(diff(win$yrange) + 1,
                                                                                  diff(win$xrange) + 1)))
       # scale such that sum{pixels} <= 1 i.e. a probability density function
       denspp           = denspp/sum(denspp)

       # calculate the probability function Fxy by normalizing denHeme to denCsr
       Fxy                  = (denspp - mean(denCsr)) / sd(denCsr)

       # cut-off based on inverse standrard normal cumulative density function
       cutoff               = qnorm(0.05, lower.tail = F)

       # with Bonferroni correction
       cutoff               = qnorm(0.05/spp$n, lower.tail = F)



       hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
       nonHotspotMask       = hotspotIm
       nonHotspotMask[nonHotspotMask > 0] = NA



       # filter out points lying outside the computed hotspot
       tmpIm                = hotspotIm
       tmpIm$v[which(tmpIm$v == 0, arr.ind = TRUE)] = NA # manually set zero pixels to NA to remove them from mask
       hotspotMask          = spatstat::as.owin(tmpIm)


       hotspotpp            = spatstat::ppp(x = spp$x,
                                            y = spp$y,
                                            window = spatstat::as.polygonal(hotspotMask),
                                            marks = spp$marks,
                                            checkdup = FALSE)


       return(list(cutoff = cutoff, # the chosen threshold above which it is a hotspot
                   GausBW = bw,		# calculated Gaussian bandwidth
                   csrspp = csr,
                   denspp = denspp,	# density image of the input point pattern
                   denCsr = denCsr,	# density image of the computed csr pattern
                   hotspotpp = hotspotpp, # the remaining points which lie within the hotspot mask
                   hotspotIm = hotspotIm, # the hotspot image
                   hotspotMask = hotspotMask, # the hotspot mask of type owin.mask
                   nonHotspotMask = nonHotspotMask)) # the non hotspot areas, to be used to hide the non hotspot areas




}


#' Calculate Gaussian bandwidth based on spatial autocorrelation - `spAutoCor`
#'
#' This function computes the Gaussian bandwidth based on sensitivity analysis of spatial autocorrlation.
#' This is used internally in `moleculaR::probMap`.
#'
#' @param spp: 	   The spatial point pattern, object of type `ppp`.
#' @param bw:        bandwidth steps at which to compute density and the corresponding Moran's I statistic.
#' @param plot:      whether to plot the result.
#' @param verbose:   whether to output progress.
#' @return
#' A numeric, the estimated Gaussian bandwidth.
#'
#' @export
#' @keywords internal

.bw.spAutoCorr <- function(spp, bw = seq(0.5,8,0.5),
                           plot = FALSE, verbose = FALSE) {


      #// create a dataframe to hold the results
      bwdf        <- data.frame(bw = bw, moransI = numeric(length(bw)))
      sppw        <- spp$marks$intensity

      if(is.null(sppw)){
            stop("spp does not have intensity weights")
      }

      #// to show progress
      if(verbose)
            pb    <- utils::txtProgressBar(min = min(bw), max = max(bw), width = 20, style = 3)



      bwdf        <- lapply(bw, function(bwi){

            if(verbose)
                  utils::setTxtProgressBar(pb, bwi)


            # create a density map for the image
            denspp      <- spatstat::density.ppp(x = spp, sigma = bwi,
                                                 weights = sppw,
                                                 W = spatstat::as.mask(spp$window, xy = list(x = spp$x, y = spp$y)))

            moranSpp    <- raster::Moran(raster::raster(denspp), w = matrix(c(1,1,1,1,0,1,1,1,1), nrow=3))


            return(data.frame(bw = bwi, moransI = moranSpp))

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
#' @param spp: 	The spatial point pattern.
#' @param weighted:  Whether the intensities are inlcude into the computation. Switch off for Collective Projections.
#' @param bw:        A vector, The gaussian band width pool used for Kernel density estimation.
#' @param csr:       Pre-computed wighted csr model corresponding to the given spp. This
#' speeds up the computation
#' @param plot:      Whether to plot.
#' @param verbose:   Whether to show progress.
#' @return
#' A list ..
#'
#' @export
#' @keywords internal
#'

.bw.iterative                 = function(spp, weighted = TRUE, bw = seq(0.5,8,0.5), csr = NULL,
                                       plot = FALSE, verbose = FALSE) {



       if(is.null(csr)) {

             ## craete a complete spatial randomness point pattern with the same number of points and window ----
             csr               = spatstat::rpoint(n = spp$n, win = spp$window)
             if(weighted){
                   csr$marks   = data.frame(intensity = sample(spp$marks$intensity))
             }
       } else { # if csr is supplied check if intensity-marks are supplied if 'weighted == TRUE'

             if(is.null(csr$marks$intensity) & weighted == TRUE)
                   stop("The supplied csr model does not contain intensity info in its marks. Consider switching 'weighted' to FALSE.\n")

       }


      if(weighted){ # if the model is intensity-weighted assign the weights
            csrw       = csr$marks$intensity
            sppw       = spp$marks$intensity
      } else { # otherwise switch weights off
            csrw       = NULL
            sppw       = NULL
      }

       #// create a dataframe to hold the results
       bwdf                   = data.frame(bw = bw, area = NA_real_)


       #// to show progress
       if(verbose)
             pb               = utils::txtProgressBar(min = min(bw), max = max(bw), width = 20, style = 3)



       bwdf$area              = sapply(bw, function(bwi){

             if(verbose)
                   utils::setTxtProgressBar(pb, bwi)



             # create a density map for csr
             denCsr               = spatstat::density.ppp(x = csr, sigma = bwi,
                                                          weights = csrw,
                                                          W = spatstat::as.mask(spp$window, xy = list(x = spp$x, y = spp$y)))

             # create a density map for the image
             denspp               = spatstat::density.ppp(x = spp, sigma = bwi,
                                                          weights = sppw,
                                                          W = spatstat::as.mask(spp$window, xy = list(x = spp$x, y = spp$y)))

             # calculate the probability function Fxy by normalizing denHeme to denCsr
             Fxy                  = (denspp - mean(denCsr)) / sd(denCsr)

             # cut-off based on inverse standrard normal cumulative density function
             cutoff               = qnorm(0.05, lower.tail = F)

             # with Bonferroni correction
             cutoff               = qnorm(0.05/spp$n, lower.tail = F)



             hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))

             hotspotIm[hotspotIm == 0] = NA # set zeros to NA to create a window
             hotspotWin          = spatstat::as.polygonal(spatstat::as.owin(hotspotIm))

             if(plot){
                   par(mfrow = c(1,1))
                   spatstat::plot.owin(spp$window, ylim = rev(spp$window$yrange), add = FALSE, main = paste0("BW = ",bwi))
                   spatstat::plot.owin(hotspotWin, col = rgb(0,1,0,1), add = TRUE)
             }


             #return(length(which(hotspotIm[,] != 0)) / prod(dim(hotspotIm)))
             return(spatstat::area.owin(hotspotWin) / spatstat::area.owin(spp$window))

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

       ep                          = .inflictPoint(x = bwdf$bw, y = bwdf$area,
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
.inflictPoint     = function(x, y, df = 7, adjRange = TRUE, plot = FALSE,
                             xQuery = seq(range(x)[1], range(x)[2], 0.1),
                             ylabel = "Moran's I") {



      # fit a smoothing spline and fine the 2nd derivative.
      smoothspl      = smooth.spline(x = x, y = y, df = df)


      deriv        = predict(smoothspl, x = xQuery, deriv = 1)$y
      spl           = predict(smoothspl, x = xQuery)$y


      if(adjRange) { # adjust the range of the  derivative to be able to plot it alongside x

            linMap <- function(i, a, b) approxfun(range(i), c(a, b))(i)

            deriv   = linMap(deriv, range(spl)[1],  range(spl)[2])

      }

      mxCurve = .getMaxCurve(xQuery, deriv)

      if(plot){
            plot(x = xQuery, y = spl, type = "b", main = "Maximum Curvature",
                 xlab = "bw",
                 ylab = ylabel)

            lines(x = xQuery, y = deriv, lty = "solid" ,  col = rgb(0,0,1,0.5), lwd = 2)
            abline(v = mxCurve, col = "chocolate1", lty ="dashed", lwd = 2)

            legend("right", bty = "n",
                   legend = c(paste0(ylabel," (bw)"), "1st-derivative",
                              paste0("infliction=", round(mxCurve, 2))
                   ),
                   col = c("black", "blue", "chocolate1"
                   ),
                   lty = c("solid", "solid","dashed"
                   ))
      }

      return(list(inflictPoint = mxCurve,
                  x = xQuery,
                  deriv = deriv,
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

.getMaxCurve                      = function(x, y, plot = FALSE)
{

      if(length(x) != length(y))
            stop("x and y have different lengths.")

      #// find the coordinates of the max bin (peak) and the last bin (tail)
      pIdx                        = which.max(y)
      tIdx                        = length(y)

      peak                        = c(x = x[pIdx], y = y[pIdx])
      tail                        = c(x = x[tIdx], y = y[tIdx])

      allPoints                   = data.frame(x = x, y = y)

      searchSpace                 = unname(unlist(apply(allPoints,
                                                        1,
                                                        FUN = function(x, p1 = peak, p2 = tail)
                                                        {
                                                              # ref: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line #
                                                              # suppose you have a line defined by two points P1 and P2, the the
                                                              # distance from point x to the line is defined (in 2D) by

                                                              abs(((p2[2] - p1[2]) * x[1]) - ((p2[1] - p1[1]) * x[2]) + (p2[1] * p1[2]) - (p2[2] * p1[1]) /
                                                                        sqrt(sum((p2 - p1) ** 2)))

                                                        })))

      searchSpace[1 : pIdx]       = NA
      maxCurveIdx                 = which.max(searchSpace)

      r                            = x[maxCurveIdx]


      if(plot)
      {
            points(rbind(peak, tail, data.frame(x = r, y = y[maxCurveIdx])),
                   col = c("red", "red", "green"), cex = 1, pch = 19)
            lines(x = rbind(peak, tail), col = "red", lty = 2)
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

       if(is.null(d)) { d <- spatstat::spatdim(X) } else check.1.integer(d)
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



