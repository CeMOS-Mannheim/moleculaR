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
                                         .spAutoCorGaussBW(spp = spp, win = win,
                                                           plot = bwPlot,
                                                           dists = dists,
                                                           numSim = numSim,
                                                           approximate = approximate,
                                                           approxRadius = approxRadius,
                                                           verbose = verbose)
                                 },
                                 iterative = {
                                         .calcGaussBW(spp = spp, win = win, weighted = weighted,
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
#' @param spp: 	   The spatial point pattern.
#' @param win:       The window object of type `spatstat::owin`.
#' @param dists:     neighborhood distances at which to compute the Moran's I statisitc..
#' @param plot:      whether to plot the result.
#' @param numSim:    the number of Monte Carlo simulations to be done, see ?spdep::moran.mc.
#' @param approximate:    when `TRUE`only points within a circle centered at the intensity-wighted
#'                         centroid (center of gravity) and radius `approxRadius` are considered. This
#'                         is used to speed up the caluclation.
#' @param approxRadius: the radius in units (usually pixels) for the approximation circlular region.By default,
#'                  this equates to the radius that produces a circular area which is half of the area of
#                  the spp window.
#' @param verbose:   whether to output progress.
#' @return
#' A numeric, the estimated Gaussian bandwidth.

#'
#' @export
#' @keywords internal

.spAutoCorGaussBW <- function(spp, win, dists = seq(1, 5, 0.5),
                              plot = FALSE, numSim = 99, approximate = FALSE,
                              approxRadius = NULL, verbose = TRUE) {

      # empty df for storage of Moran's I values
      moranI <- data.frame(d = numeric(0), moranStat = numeric(0), pvalue = numeric(0))

      # to speed things up
      if(approximate) {

            # define focus area
            centroid <- c(weighted.mean(spp$x, spp$marks$intensity), c(weighted.mean(spp$y, spp$marks$intensity)))
            approxRadius <- ifelse(is.null(approxRadius), sqrt(spatstat::area(win)/(2*pi)), approxRadius)
            focusWin <- spatstat::disc(radius = approxRadius, centroid)


            focusppp            = spatstat::ppp(x = spp$x,
                                                y = spp$y,
                                                window = focusWin,
                                                marks = spp$marks,
                                                checkdup = FALSE)

            focusppp  <- spatstat::as.ppp(focusppp)

            # coordinates and intensities
            xyCoords <- as.matrix(spatstat::coords(focusppp))
            intensities <- focusppp$marks$intensity

      } else {
            # coordinates and intensities
            xyCoords <- as.matrix(spatstat::coords(spp))
            intensities <- spp$marks$intensity

      }



      #// to show progress
      if(verbose)
            pb    <- utils::txtProgressBar(min = min(dists), max = max(dists), width = 20, style = 3)

      # loop d through the distances

      moranI      <- lapply(dists, function(di){
            if(verbose)
                  utils::setTxtProgressBar(pb, di)

            nb    <- spdep::dnearneigh(xyCoords, d1 = 0, d2 = di)
            lw    <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
            mcs   <- spdep::moran.mc(intensities, lw, nsim = numSim, zero.policy = TRUE)

            data.frame(d = di, moranStat = mcs$statistic, pvalue = mcs$p.value)

      })

      moranI      <- do.call("rbind", moranI)

      if(plot){
            plot(x = moranI$d, y = moranI$moranStat, type = "b", main = "Sensitivity analsis",
                 xlab = "Neigghborhood distances (pixels)",
                 ylab = "Moran's I statistic")

            abline(v =  moranI$d[which.max(moranI$moranStat)], col = "chocolate1", lty ="dashed", lwd = 2)

            if(approximate){

                  spatstat::plot.owin(win, ylim = rev(win$yrange))
                  spatstat::plot.ppp(focusppp, which.marks = "intensity", cols = viridis::viridis_pal(option = "inferno")(10), add = T)
                  plot(focusWin, add = TRUE, do.col=T)

            }


      }

      if(verbose)
            close(pb); cat("Done. \n")

     #return(list(moranIStats = moranI, bw = moranI$d[which.max(moranI$moranStat)]))
      return(c(bw = moranI$d[which.max(moranI$moranStat)]))

}


#' Calculate Gaussian Bandwidth - `iterative`
#'
#' This function calculates the Gaussian bandwidth to be used for analyte probability maps by
#' finding the curve infliction point for a curve that represents hotspot area as a function
#' of bandwidth. This is used internally in `moleculaR::probMap`.
#'
#' @param spp: 	The spatial point pattern.
#' @param win:       The window object of type `spatstat::owin`.
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

.calcGaussBW                 = function(spp, win, weighted = TRUE, bw = seq(1,5,0.5), csr = NULL,
                                       plot = TRUE, verbose = TRUE) {



       if(is.null(csr)) {

             ## craete a complete spatial randomness point pattern with the same number of points and window ----
             csr               = spatstat::rpoint(n = spp$n, win = win)
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
                                                          W = spatstat::as.mask(win,
                                                                                dimyx=c(diff(win$yrange) + 1,
                                                                                        diff(win$xrange) + 1)))

             # create a density map for the image
             denspp               = spatstat::density.ppp(x = spp, sigma = bwi,
                                                          weights = sppw,
                                                          W = spatstat::as.mask(win,
                                                                                dimyx=c(diff(win$yrange) + 1,
                                                                                        diff(win$xrange) + 1)))

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
                   spatstat::plot.owin(win, ylim = rev(win$yrange), add = FALSE, main = paste0("BW = ",bwi))
                   spatstat::plot.owin(hotspotWin, col = rgb(0,1,0,1), add = TRUE)
             }


             #return(length(which(hotspotIm[,] != 0)) / prod(dim(hotspotIm)))
             return(spatstat::area.owin(hotspotWin) / spatstat::area.owin(win))

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

       ep                          = .inflictPoint(x = bwdf$bw, y = bwdf$area)


      if(plot) {

             plot(x = bwdf$bw, y = bwdf$area, type = "p", main = "MPM significance area to total tissue area",
                  xlab = "Gaussian Bandwidth",
                  ylab = "MPM area / total area")

             lines(x = ep$x, y = ep$spl, lty = "solid" ,  col = rgb(0,0,0,0.5), lwd = 2)
             lines(x = ep$x, y = ep$deriv2,lty = "solid" ,  col = rgb(0,0,1,0.5), lwd = 2)
             abline(v =  ep$inflictPoint, col = "chocolate1", lty ="dashed", lwd = 2)

             legend("bottomright",
                    legend = c("area(bw)", "2nd-derivative",
                               paste0("optimum bw=", round(ep$inflictPoint, 2))#,
                               #paste0("ppl bw=", round(bwppl, 2)),
                               #paste0("scott bw=", round(bwscott, 2))
                               ),
                    col = c("black", "blue", "chocolate1"#,
                            #"green",
                            #"red"
                            ),
                    lty = c("solid", "solid","dashed"#,
                            #"dashed", "dashed"
                            ))

      }

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
#' @param xQuery:    x values, at which to evaluate the fitted spline.
#' A list ..
#'
#' @export
#' @keywords internal
#'
.inflictPoint     = function(x, y, df = 7, adjRange = TRUE,
                        xQuery = seq(range(x)[1], range(x)[2], 0.1)) {



       # fit a smoothing spline and fine the 2nd derivative.
       smoothspl      = smooth.spline(x = x, y = y, df = df)


       deriv2        = predict(smoothspl, x = xQuery, deriv = 2)$y
       spl           = predict(smoothspl, x = xQuery)$y


       if(adjRange) { # adjust the range of the 2nd derivative to be able to plot it alongside x

              linMap <- function(x, a, b) approxfun(range(x), c(a, b))(x)

              deriv2   = linMap(deriv2, range(y)[1],  range(y)[2])

       }

       return(list(inflictPoint = xQuery[which.min(deriv2)],
                   x = xQuery,
                   deriv2 = deriv2,
                   spl = spl))




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



