#' Calculate molecular probablity maps
#'
#' This function creates molecular probability maps from a given spatial point pattern representation.
#'
#' @param sppMoi 	   The spatial point pattern.
#' @param bw        The Gaussian bandwidth used for Kernel density estimation. If a single numeric value is supplied
#' then it will be used for the computation of KDE. If a vector is supplied (in increasing order), then these values
#' will be passed to internal method `spAutoCor` for computing the bandwidth based on cross validation,
#' see `?moleculaR::.bw.spAutoCorr`.
#' @param edgeCorrection   a logical. Whether to apply edge correction through
#' erosion by calling `spatstat.geom::erosion`. This could be beneficial for
#' removing edge artifacts which might otherwise impact the result.
#' @param bwMethod The method used for computing the Gaussian bandwidth, c("spAutoCor", "scott").
#' @param diagPlots    whether to plot the Gaussian bw selection procedure, ignored when \code{bwMethod = "scott"}.
#' @param csrIntensities How the intensities for the csrMoi model are generated, c("resample", "Poisson", "Gaussian").
#' @param reference   An sppMoi object that is designated as the "control" or "null" alternative of \code{sppMoi}. If supplied the
#' csrMoi model will be generated from the intensities of this object.
#' @param csrMoi       Pre-computed wighted csrMoi model corresponding to the given sppMoi. This
#' could be useful when generating hotspot-bw curves. For internal use only.
#' @param pvalThreshold The p-value threshold to be used for the hypothesis testing.
#' @param pvalCorrection The method used for p-values correction, see '?p.adjust' for more details.
#' @param seed    a single value, interpreted as integer or `NULL` (default) controlling the process of random
#' number generation. If set, ensures that the same CSR model is generated at different runs for the
#' same input `sppMoi`
#' @param verbose   whether to output progress.
#' @param ...: arguments passed to `plot` for bandwidth plotting when `bwMethod="spAutoCor"`.

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
#' @include manualSpatstatImport.R
#'
probMap                     <- function(sppMoi,
                                        reference = NULL,
                                        bw = seq(0.5, 10),
                                        edgeCorrection = FALSE,
                                        bwMethod = "spAutoCor",
                                        csrMoi = NULL,
                                        pvalThreshold = 0.05,
                                        pvalCorrection = "BH",
                                        diagPlots = FALSE,
                                        seed = NULL,
                                        verbose = TRUE,
                                        ...) {



      if(!("analytePointPattern" %in% class(sppMoi))){
            stop("sppMoi must be of 'analytePointPattern' and 'ppp' class. \n")
      }

      if(sppMoi$n < 100){
            warning("The supplied sppMoi has too few points (detected peaks). ",
                    "Be advised that this might negatively affect the subsequent spatial analysis/results. \n")
      }

      if(is.null(reference)){

            is.crossTissueCase <- FALSE

      } else {

            is.crossTissueCase <- TRUE

            if(!(is.ppp(reference))) {
                  stop("Supplied reference is not an ppp object. \n")
            }

            if(!("analytePointPattern" %in% class(reference))){
                  stop("reference must be of 'analytePointPattern' and 'ppp' class. \n")
            }

            if(is.null(reference$marks$intensity)){
                  stop("The supplied reference does not contain intensity info in its marks.\n")
            }

            if(reference$n == 0){
                  stop("The supplied reference is empty (n = 0).\n")
            }

            if(reference$n < 100){
                  warning("The supplied reference has too few points (detected peaks). ",
                          "Be advised that this might negatively affect the subsequent eCDF analysis/results. \n")
            }

      }

      set.seed(seed)

      if(edgeCorrection){ # apply slight erosion to edges

            sppMoi <- .erosion(sppMoi)

      }

      if(is.null(csrMoi)) { # if a csr model was not provided

            csrMoi <- .createCSR(sppMoi)

      } else { # if csrMoi is supplied check if intensity-marks are supplied

            if(!("analytePointPattern" %in% class(csrMoi))){
                  stop("csrMoi must be of 'analytePointPattern' and 'ppp' class. \n")
            }

            if(is.null(csrMoi$marks$intensity)){
                  stop("The supplied csrMoi model does not contain intensity info in its marks.\n")
            }


      }


      # compute kernel density bandwidth
      if(length(bw) > 1)
      {

            if(any(diff(bw) < 0)){
                  stop("supplied bw must be either a single numeric of a vector of numerics of increasing order.\n")
            }


            bw        <- switch(bwMethod,
                                spAutoCor = {
                                      .bw.spAutoCorr(sppMoi = sppMoi,
                                                     plot = diagPlots,
                                                     bw = bw,
                                                     verbose = verbose,
                                                     ...)
                                },
                                scott = {
                                      .bw.scott.iso2(sppMoi)
                                },
                                stop("wrong bwMethd specified, must be one of c('spAutoCor','scott')"))



      }



      ## __ statistical testing and pvalue correction __ #

      # probability density function - csr
      rhoCsr        <- .rho(spp = csrMoi, sigma = bw)
      # probability density function - MOI
      rhoMoi        <- .rho(spp = sppMoi, sigma = bw)

      ## hypthothesis testing for spatial autocorrelation ##

      # probability values for coldspots (lower tail test) - coldspot
      pvalsLwrSpatial <- .spatialHTest(rhoCsr, rhoMoi, pvalCorrection, TRUE)
      # probability values for hotspot (upper tail test) - hotspot
      pvalsUprSpatial <- .spatialHTest(rhoCsr, rhoMoi, pvalCorrection, FALSE)



      # if a reference spp was provided an additional test should be carried out
      # to test the liklihood of each intensity in the test tissue belonging to
      # the distribution of the control (reference) tissue intensities.
      if(is.crossTissueCase){


            ## hypthothesis testing for intensities of test vs reference intensities ##

            # lower tail - coldspot
            pvalsLwrIntensity <- .intensHTest(reference, sppMoi, pvalCorrection, TRUE)
            # upper tail - hotspot
            pvalsUprIntensity <- .intensHTest(reference, sppMoi, pvalCorrection, FALSE)

            if(diagPlots){

                  print(.ecdfPlot(reference, sppMoi))
                  print(.bvPlot(reference, sppMoi))
            }

      }


      ## __ hotspot __ ##

      hotspotIm        <- eval.im(rhoMoi * (pvalsUprSpatial <= pvalThreshold))

      if(is.crossTissueCase){ # points of sppMoi that are not drawn from the distribution of reference
            hotspotIm        <- eval.im(hotspotIm * (pvalsUprIntensity <= pvalThreshold))
      }

      # filter out points lying outside the computed hotspot
      tmpIm            <- hotspotIm
      tmpIm$v[which(tmpIm$v == 0, arr.ind = TRUE)] <- NA # manually set zero pixels to NA to remove them from mask
      hotspotMask      <- as.owin(tmpIm)


      hotspotpp        <- ppp(x = sppMoi$x,
                              y = sppMoi$y,
                              window = as.polygonal(hotspotMask),
                              marks = sppMoi$marks,
                              checkdup = FALSE)

      hotspotpp        <- analytePointPattern(spp = hotspotpp, mzVals = sppMoi$metaData$mzVals,
                                              metaData = as.list(sppMoi$metaData)) # only for conformity with S3 class


      ## __ coldspot __ ##

      coldspotIm        <- eval.im(rhoMoi * (pvalsLwrSpatial <= pvalThreshold))

      if(is.crossTissueCase){ # points of sppMoi that are not drawn from the distribution of reference
            coldspotIm        <- eval.im(coldspotIm * (pvalsLwrIntensity <= pvalThreshold))
      }

      # filter out points lying outside the computed coldspot
      tmpIm             <- coldspotIm
      tmpIm$v[which(tmpIm$v == 0, arr.ind = TRUE)] <- NA # manually set zero pixels to NA to remove them from mask
      coldspotMask      <- as.owin(tmpIm)


      coldspotpp        <- ppp(x = sppMoi$x,
                               y = sppMoi$y,
                               window = as.polygonal(coldspotMask),
                               marks = sppMoi$marks,
                               checkdup = FALSE)

      coldspotpp        <- analytePointPattern(spp = coldspotpp, mzVals = sppMoi$metaData$mzVals,
                                               metaData = as.list(sppMoi$metaData)) # only for conformity with S3 class


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

.createCSR <- function(sppMoi, reference = NULL) {

      # sppMoi: The spatial point pattern under investigation.
      # reference: The spatial point pattern of a control (reference) tissue. -> Depricated

      if(is.null(reference)){

            csrMoi         <- rpoint(n = sppMoi$n, win = sppMoi$window)

            csrMoi$marks   <- data.frame(intensity = sample(sppMoi$marks$intensity))

            csrMoi <- analytePointPattern(spp = csrMoi, mzVals = sppMoi$metaData$mzVals,
                                          metaData = as.list(sppMoi$metaData))

            return(csrMoi)

      }else{

            csrMoi              <- rpoint(n = sppMoi$n, win = sppMoi$window)


            if(reference$n >= csrMoi$n){
                  csrMoi$marks <- data.frame(intensity = sample(reference$marks$intensity,
                                                                size = csrMoi$n,
                                                                replace = FALSE))
            }else{
                  csrMoi$marks <- data.frame(intensity = sample(reference$marks$intensity,
                                                                size = csrMoi$n,
                                                                replace = TRUE))
            }


            csrMoi <- analytePointPattern(spp = csrMoi, mzVals = sppMoi$metaData$mzVals,
                                          metaData = as.list(sppMoi$metaData))
            return(csrMoi)

      }

}

# computes the spatial density function rho
.rho <- function(spp, sigma){
      # spp: a spatial point pattern
      # sigma: bw for the KDE

      # extract weights - intensities
      w       <- spp$marks$intensity

      # create a spaital window
      win <- as.mask(spp$window,
                     dimyx=c(diff(spp$window$yrange) + 1,
                             diff(spp$window$xrange) + 1))

      # create a density map for spp
      sppRho  <- density.ppp(x = spp, sigma = sigma,
                             weights = w,
                             W = win, positive = TRUE)


      # scale such that sum{pixels} <= 1 i.e. a probability density function
      sppRho           <- sppRho/sum(sppRho)

      return(sppRho)


}


.spatialHTest <- function(rhoCsr, rhoMoi, pvalCorrection, lower.tail){
      # Performs hypothesis testing for each intensity in rhoMoi (probability
      # of it being drawn from the distribution of rhoCsr)
      #
      # rhoCSR: the spatial density function of the csr spp
      # rhoMOI: the spatial density function of the MOI spp
      # pvalCorrection: the p-value correction method as in p.adjust.
      # lower.tail: logical, lower or upper tail?
      #
      # returns an image whose intensities are probabilities.

      # null hypothesis
      mucsrMoi            <- mean(rhoCsr, na.rm = TRUE)
      sigmacsrMoi         <- sd(rhoCsr, na.rm = TRUE)

      # convert to data.frame
      rhoMoidf <- as.data.frame.im(x = rhoMoi)
      pvals <- rhoMoidf

      # generate p-values - lower tail or upper tail
      pvals$value <- pnorm(rhoMoidf$value, mean = mucsrMoi, sd = sigmacsrMoi, lower.tail = lower.tail)

      # correction
      pvals$value <- p.adjust(p = pvals$value, method = pvalCorrection)

      # convert back to image
      pvals <- as.im.data.frame(pvals)

      return(pvals)

}


.intensHTest <- function(sppRef, sppTest, pvalCorrection, lower.tail){
      # Performs hypothesis testing for each intensity in imTest (probability
      # of it being drawn from the distribution of imRef) based on eCDF. This
      # only performed when control (refernce) tissue is available.
      #
      # sppRef: the image created of the reference tissue sppMoi.
      # sppTest: the image created of the test tissue sppMoi.
      # pvalCorrection: the p-value correction method as in p.adjust.
      # lower.tail: logical, lower or upper tail?
      #
      # returns an image whose intensities are probabilities.

      # convert to data.frame
      dfRef <- as.data.frame.im(spp2im(sppRef, rescale = FALSE, zero.rm = TRUE))
      dfTest <- as.data.frame.im(spp2im(sppTest, rescale = FALSE, zero.rm = TRUE))
      pvals <- dfTest


      # null hypothesis - reference (control) tissue
      ecdfRef <- ecdf(dfRef$value)

      if(lower.tail){

            # lower tail - coldspot
            pvals$value <- ecdfRef(dfTest$value)
            pvals$value <- p.adjust(pvals$value, method = pvalCorrection)


      } else {

            # upper tail - hotspot
            pvals$value <- 1-ecdfRef(dfTest$value)
            pvals$value <- p.adjust(pvals$value, method = pvalCorrection)

      }



      # add a small number to zero probability to distinguish then from empty pixels
      # zridx <- which(pvals$value == 0)
      # pvals$value[zridx] <- 1e-10

      # convert to image
      pvals <- as.im.data.frame(pvals)

      # set empty pixels to NA
      #pvals$v[which(pvals$v == 0, arr.ind = TRUE)] <- NA

      return(pvals)

}

.bvPlot <- function(sppRef, sppTest){

      plotlist <- list(Test = sppTest$marks$intensity,
                       Reference = sppRef$marks$intensity)

      r <- range(unlist(plotlist), na.rm = TRUE)
      plotmx <- r[2] + 4*(diff(r)/10)
      txtmx <- r[2] + 3*(diff(r)/10)
      segmenentmx <- r[2] + 2*(diff(r)/10)
      tickmx <- r[2] + 1*(diff(r)/10)
      mn <- r[1]

      par(bty = "n",  cex.axis = 1.5, cex.lab = 1.5, mar = c(5.1, 5.1, 4.1, 2.1))
      vioplot::vioplot(x = plotlist, ylab ="Intensities", main = "Cross-tissue Intensities",
                       col=c(to.transparent("#FF8C69", 0.5), to.transparent("#00FFFF", 0.5)),
                       plotCentre = "point", rectCol = "white", colMed = "black",
                       areaEqual = F, wex = 1, ylim = c(mn,plotmx))
      # line segments
      segments(x0 = 1, y0 = segmenentmx, x1 = 2, y1 = segmenentmx, col = "black")
      segments(x0 = 1, y0 = segmenentmx, x1 = 1, y1 = tickmx, col = "black")
      segments(x0 = 2, y0 = segmenentmx, x1 = 2, y1 = tickmx, col = "black")


      #text(x=1.5,y=max(unlist(plotlist))+0.1,"***",pos=3,cex=1.5)

      wt <- wilcox.test(x = plotlist$Test, y = plotlist$Reference)

      text(x=1.5,y=txtmx,"Mann-Whitney U Test", pos=3, cex=1)

      if(wt$p.value < 0.001){
            text(x=1.5,y=segmenentmx,"p-value < 0.001", pos=3, cex=1)
      } else {
            text(x=1.5,y=segmenentmx,paste0("p-value = ", round(wt$p.value, 2)), pos=3, cex=1)
      }

}


.ecdfPlot <- function(sppRef, sppTest){

      plotlist <- list(Test = sppTest$marks$intensity,
                       Reference = sppRef$marks$intensity)

      ecdfTest <- ecdf(plotlist$Test)
      ecdfRef <- ecdf(plotlist$Reference)
      par(bty = "n",  cex.axis = 1.5, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1))
      plot(ecdfTest, xlim = range(unlist(plotlist), na.rm = TRUE),
           ylab = "eCDF", xlab = "Intensity Quantiles", main = "Cross-tissue eCDFs",
           bty = "n", col = "#FF8C69", lwd = "2")
      lines(ecdfRef, col = "#00FFFF", lwd = "2", lty = "dashed")
      legend("topleft", legend = c("Test", "Reference"), col = c("#FF8C69", "#00FFFF"),
             lty = c("solid", "dashed"), bty ="n")


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
#' @param ...: arguments passed to `plot` for bandwidth plotting.
#' @return
#' A numeric, the estimated Gaussian bandwidth.
#'
#' @export
#' @keywords internal

.bw.spAutoCorr <- function(sppMoi, bw = c(0.5, seq(1, 9, 1)),
                           plot = FALSE, verbose = FALSE, ...) {


      #// create a dataframe to hold the results
      bwdf        <- data.frame(bw = bw, moransI = numeric(length(bw)))
      win         <- as.mask(sppMoi$window, dimyx=c(diff(sppMoi$window$yrange) + 1,
                                                    diff(sppMoi$window$xrange) + 1))


      if(is.null(sppMoi$marks$intensity)){
            stop("sppMoi does not have intensity weights")
      }

      #// to show progress
      if(verbose)
            pb    <- utils::txtProgressBar(min = min(bw), max = max(bw), width = 20, style = 3)



      bwdf        <- lapply(bw, function(bwi){

            if(verbose)
                  utils::setTxtProgressBar(pb, bwi)




            # create a density map for the image
            rhoMoi      <- density.ppp(x = sppMoi, sigma = bwi,
                                       weights = sppMoi$marks$intensity,
                                       W = win)



            moransppMoi    <- raster::Moran(raster::raster(rhoMoi))


            return(data.frame(bw = bwi, moransI = moransppMoi))

      })

      if(verbose)
            close(pb)

      bwdf <- do.call("rbind", bwdf)


      #// compute the elbow point
      ep                <- kneePoint(x = bwdf$bw, y = bwdf$moransI,
                                     plot = plot, ...)


      if(verbose)
            cat("done.\n")

      return(c(bw = ep))

}




#' Find Knee (or elbow) point of a curve
#'
#' This function calculates the knee/elbow point of a curve based on the kneedle
#' algorithm (satopaa et al, 2011). This is used internally in `moleculaR::.calcGaussBW`.
#' This is a simplified implementation.
#'
#' @param x: 	   x values representing the bandwidth values
#' @param y:         y values representing the Moran's I statistic.
#' @param df:        degrees of freedom for the smoothing spline.
#' @param plot:   whether to plot the result.
#' @param xQuery:    x values to be used for smoothing the original curve via
#' a fitted spline.
#' @param sign +1 for increasing values (knee) and  -1 for decreasing values (elbow).
#' @param ...: arguments passed to `plot`.
#'
#' @return
#' Returns a numeric, the calculated knee point representing the optimum bandwidth.
#'
#' @export
#' @keywords internal
#'
#' @references
#' Satopaa, Ville, et al. "Finding a" kneedle" in a haystack: Detecting knee points
#' in system behavior." 2011 31st international conference on distributed computing
#' systems workshops. IEEE, 2011. (doi: 10.1109/ICDCSW.2011.20)
#'


kneePoint     <- function(x, y, df = 7,
                          xQuery = seq(range(x)[1], range(x)[2], 0.1),
                          plot = FALSE, sign = +1, ...) {




      # fit a smoothing spline/loess
      smoothx   <- xQuery
      smoothspl <- smooth.spline(x = x, y = y, df = df)
      smoothy <- predict(smoothspl, x = smoothx)$y


      # normalize points of the smoothing curve to unit square
      smoothnx     <- (smoothx - min(smoothx)) / (max(smoothx) - min(smoothx))
      smoothny     <- (smoothy - min(smoothy)) / (max(smoothy) - min(smoothy))


      # apply kneedle
      k <- .kneedle(x = smoothnx, y = smoothny, sign = sign)


      if(plot){
            par(mar=c(5.1, 5.1, 4.1, 2.1))
            plot(x = smoothx, y = smoothy, type = "l",
                 main = "Knee-point estimation",
                 xlab = "Gaussian Bandwidth",
                 ylab = "Moran's I",
                 ...)

            # defaults
            lwd <- 1; cex <- 1

            # extract size arg from ...
            pargs <- list(...)

            if(length(pargs) > 0){
                  n <- names(pargs)

                  if("lwd" %in% n){
                        lwd <- pargs$lwd
                  }


                  if("cex.lab" %in% n){
                        cex <- pargs$cex.lab
                  }
            }

            abline(v = smoothx[k], col = "chocolate1", lty ="dashed", lwd = lwd)

            legend("right", bty = "n",
                   legend = c(paste0("Moran's I"),
                              paste0("Knee pnt=", round(smoothx[k], 2))
                   ),
                   col = c("black",
                           "chocolate1"
                   ),
                   lty = c("solid",
                           "dashed"
                   ),
                   lwd = lwd, cex = cex)
      }

      return(smoothx[k])


}


#' Find the knee/elbow point in a vector using the Kneedle algorithm.
#'
#' This function uses the Kneedle algorithm to find the index of the knee point
#' in the provided x,y-vectors.
#'
#' @param x numeric vector, x-values.
#' @param y numeric vector, y-values.
#' @param sign +1 for increasing values and  -1 for decreasing values.
#'
#' @return The index of the knee/elbow.
#'
.kneedle <- function(x, y, sign = 1) {

      if(length(x) != length(y))
            stop("error in internal function .kneedle; x and y of different lengths.\n")

      start = c(x[1], y[1])
      end = c(x[length(x)], y[length(y)])

      k <- which.max(lapply(1:length(x), function(i) {
            sign * -1 * .dist2d(c(x[i], y[i]),
                                start,
                                end)
      }))


      k

}


.dist2d <- function(a,b,c) {
      v1 <- b - c
      v2 <- a - b
      m <- cbind(v1,v2)
      d <- det(m)/sqrt(sum(v1*v1))
      d
}


#'  Maximum curvature
#'
#' Finds the maximum curvature (elbow point) of vector of values representing a curve by finding the maximum distance from the
#' curve to the line drawn between the peak and tail of that curve. **Deprecated**.
#'
#' @param x:    a numeric vector.
#' @param plot: a logical to plot the result for validation purposes. Defaults to FALSE.
#'
#'
#' @return Returns the maximum curvature threshold (elbow point) of `x`.
#'
#'
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
#' Note that this function is already available in recent \code{spatstat} versions (see \code{?.bw.scott}).
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

      if(is.null(d)) { d <- spatdim(X) } else check.1.integer(d)
      nX <- npoints(X)
      cX <- coords(X, spatial=TRUE, temporal=FALSE, local=FALSE)
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



#' Apply erosion to an spp
#'
#' This function is equivalent to `spatstat.geom::erosion` but applied to the `spp`
#' object and retains its attributes
#'
#' @param spp: 	   A spatial point pattern (object of class \code{spp}).
#' @return
#' Another `spp` object with its window eroded on the edges with a fixed r = 0.5
#' length units (pixels).
#'
.erosion <- function(spp){

      if (!("analytePointPattern" %in% class(spp))) {
            stop(".erosion: spp must be of 'analytePointPattern' and 'ppp' class. \n")
      }
      erodedWin <- erosion(w = spp$window, r = 0.5)
      spp <- analytePointPattern(x = spp$x, y = spp$y, intensity = spp$marks$intensity,
                                 win = erodedWin, mzVals = spp$metaData$mzVals, metaData = as.list(spp$metaData))
      return(spp)

}


