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
#' @return
#' A list ..
#'
#' @export
#'
calcGaussBW                 = function(spp, win, weighted = TRUE, bw = seq(1,12,0.5), csr = NULL, plot = TRUE) {



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

       #// create a dataframe to hold the results
       bwdf                   = data.frame(bw = bw, area = NA_real_)

       cat("\n")
       for(i in 1 : length(bw)) {

              cat("\r", "i = ", i, "of", length(bw), "              ")

             if(weighted){ # if the model is intensity-weighted assign the weights
                   csrw       = csr$marks$intensity
                   sppw       = spp$marks$intensity
             } else { # otherwise switch weights off
                   csrw       = NULL
                   sppw       = NULL
             }

              # create a density map for csr
              denCsr               = spatstat::density.ppp(x = csr, sigma = bw[i],
                                                           weights = csrw,
                                                           W = spatstat::as.mask(win,
                                                                                 dimyx=c(diff(win$yrange) + 1,
                                                                                         diff(win$xrange) + 1)))

              # create a density map for the image
              denspp               = spatstat::density.ppp(x = spp, sigma = bw[i],
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
              #nonHotspotMask       = hotspotIm
              #nonHotspotMask[nonHotspotMask > 0] = NA


              bwdf$area[i]         = length(which(hotspotIm[,] != 0)) / prod(dim(hotspotIm))

       }
       cat("\n")

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

             plot(x = bwdf$bw, y = bwdf$area, type = "p", main = "Hotspot area to total tissue area(bw)", xlab = "bw",
                  ylab = "Hotspot to total tissue area")

             lines(x = ep$x, y = ep$spl, lty = "solid" ,  col = rgb(0,0,0,0.5))
             lines(x = ep$x, y = ep$deriv2,lty = "solid" ,  col = rgb(0,0,1,0.5))
             abline(v =  ep$elbowPoint, col = "gold", lty ="dashed")
             #abline(v =  bwppl, col = "green", lty ="dashed")
             #abline(v =  bwscott, col = "red", lty ="dashed")
             #spatstat::bw.ppl = 3.981072
             #spatstat::bw.scott = 8.244953
             legend("bottomright",
                    legend = c("area(bw)", "2nd-derivative",
                               paste0("spline bw=", round(ep$elbowPoint, 2))#,
                               #paste0("ppl bw=", round(bwppl, 2)),
                               #paste0("scott bw=", round(bwscott, 2))
                               ),
                    col = c("black", "blue", "gold"#,
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


