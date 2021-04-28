#' Calculate Gaussian Bandwidth
#'
#' This function calculates the Gaussian bandwidth to be used for analyte probability maps by
#' finding the curve infliction point for a curve that represents hotspot area as a function
#' of bandwidth.
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
calcGaussBW                 = function(spp, win, weighted = TRUE, bw = seq(1,10,0.5), csr = NULL,
                                       plot = TRUE, verbose = TRUE) {



       if(is.null(csr)) {

             ## craete a complete spatial randomness point pattern with the same number of points and window ----
             csr               = spatstat::rpoispp(lambda = spatstat::intensity(spp), win = win)
             if(weighted){
                   csr$marks   = data.frame(intensity = rpois(csr$n, median(spp$marks$intensity)))
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

       # for(i in 1 : length(bw)) {
       #
       #       utils::setTxtProgressBar(pb, i)
       #
       #
       #
       #        # create a density map for csr
       #        denCsr               = spatstat::density.ppp(x = csr, sigma = bw[i],
       #                                                     weights = csrw,
       #                                                     W = spatstat::as.mask(win,
       #                                                                           dimyx=c(diff(win$yrange) + 1,
       #                                                                                   diff(win$xrange) + 1)))
       #
       #        # create a density map for the image
       #        denspp               = spatstat::density.ppp(x = spp, sigma = bw[i],
       #                                                     weights = sppw,
       #                                                     W = spatstat::as.mask(win,
       #                                                                           dimyx=c(diff(win$yrange) + 1,
       #                                                                                   diff(win$xrange) + 1)))
       #
       #        # calculate the probability function Fxy by normalizing denHeme to denCsr
       #        Fxy                  = (denspp - mean(denCsr)) / sd(denCsr)
       #
       #        # cut-off based on inverse standrard normal cumulative density function
       #        cutoff               = qnorm(0.05, lower.tail = F)
       #
       #        # with Bonferroni correction
       #        cutoff               = qnorm(0.05/spp$n, lower.tail = F)
       #
       #
       #
       #        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
       #        #nonHotspotMask       = hotspotIm
       #        #nonHotspotMask[nonHotspotMask > 0] = NA
       #
       #
       #        bwdf$area[i]         = length(which(hotspotIm[,] != 0)) / prod(dim(hotspotIm))
       #
       # }
       # close(pb)
       # cat("\n")

       #// compute the elbow point
       if(all(bwdf$area == 0)) { # no hotspot detected

              cat("no hotspot detected for all bw. Setting arbitrary bw .. \n")
              return(list(bwdf = bwdf,
                          elbowPointData = list(elbowPoint = 3))) # arbitrary bw

       }

       ep                          = elbowPointSimple(x = bwdf$bw, y = bwdf$area)
       #bwppl                       = spatstat::bw.ppl(spp)
       #bwscott                     = bw.scott.iso(spp)

      if(plot) {

             plot(x = bwdf$bw, y = bwdf$area, type = "p", main = "MPM significance area to total tissue area",
                  xlab = "Gaussian Bandwidth",
                  ylab = "MPM area / total area")

             lines(x = ep$x, y = ep$spl, lty = "solid" ,  col = rgb(0,0,0,0.5), lwd = 2)
             lines(x = ep$x, y = ep$deriv2,lty = "solid" ,  col = rgb(0,0,1,0.5), lwd = 2)
             abline(v =  ep$elbowPoint, col = "chocolate1", lty ="dashed", lwd = 2)
             #abline(v =  bwppl, col = "green", lty ="dashed")
             #abline(v =  bwscott, col = "red", lty ="dashed")
             #spatstat::bw.ppl = 3.981072
             #spatstat::bw.scott = 8.244953
             legend("bottomright",
                    legend = c("area(bw)", "2nd-derivative",
                               paste0("optimum bw=", round(ep$elbowPoint, 2))#,
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

       return(list(bwdf = bwdf,
                   elbowPointData = ep)) # the non hotspot areas, to be used to hide the non hotspot areas




}


