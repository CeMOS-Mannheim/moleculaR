#' Calculate molecular probablity maps
#'
#' This function creates molecular probability maps from a given spatial point pattern representation.
#'
#' @param spp: 	The spatial point pattern.
#' @param win:       The window object of type `spatstat::owin`.
#' @param weighted:  Whether the intensities are inlcude into the computation.
#' @param bw:        The Gaussian bandwidth used for Kernel density estimation.
#' @param bwMethod: The method used for computing the Gaussian bandwidth, c("iterative", "scott").
#' @param bwPlot:    whether to plot the Gaussian bw selection procedure, ignored when \code{bwMethod = "scott"}.
#' @param csrIntensities: How the intensities for the CSR model are generated, c("Poisson", "Gaussian"). 
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
probMap                     = function(spp, win, weighted = TRUE, bw = NULL, control = NULL,
                                       bwMethod = "iterative", csrIntensities = "Poisson", csr = NULL, bwPlot = FALSE,
                                       verbose = TRUE) {




       if(is.null(csr)) {

              ## craete a complete spatial randomness point pattern with the same number of points and window ----

             if(is.null(control)){ # if control is not supplied generate the csr model from spp

                   csr               = spatstat::rpoispp(lambda = spatstat::intensity(spp), win = win)
                   if(weighted){
                         csr$marks   = switch(csrIntensities,
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

                   if(weighted){
                         if(is.null(control$marks$intensity) & weighted == TRUE){
                               stop("The supplied control does not contain intensity info in its marks.\n")
                               }
                         csr$marks   = switch(csrIntensities,
                                                "Poisson" = {
                                                      data.frame(intensity = rpois(csr$n, mean(control$marks$intensity)))
                                                },
                                                "Gaussian"= {
                                                      data.frame(intensity = rnorm(csr$n, mean(control$marks$intensity),
                                                                                    sd(control$marks$intensity)))
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
                                 iterative = {
                                         calcGaussBW(spp = spp, win = win, weighted = weighted,
                                                     csr = csr, plot = bwPlot,
                                                     verbose = verbose)$elbowPointData$elbowPoint
                                 },
                                 scott = {
                                         .bw.scott.iso(spp)
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

       # create a density map for the image
       denspp               = spatstat::density.ppp(x = spp, sigma = bw,
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


.bw.scott <- function(X, isotropic=FALSE, d=NULL) {
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

.bw.scott.iso <- function(X) { .bw.scott(X, isotropic=TRUE) }

