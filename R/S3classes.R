## Classes
##
## Under the hood, moleculaR relies on MALDIquant, MALDIquantForeign,
## Matrix and spatstat packages.It inherits classes from 'Matrix' and
## 'spatstat' packages. For more info please check the documentation
## the respective packages.
##
##
## to document this -> https://r-pkgs.org/man.html#man-functions

#' @include manualSpatstatImport.R
NULL

#' fwhm Class Constructor
#'
#' This constructs a 'fwhm' object which stores information about the calculated full width half maximum.
#' This function is for internal use and normally should not be called by the user.
#'
#' @param fwhmInterpolator A function that is used to find fwhm value for any given m/z.
#' @param  peaks  The m/z values of the peaks of a single spectrum.
#' @param  fwhmVals   The corresponding fwhm values in Da.
#'
#' @return
#' An S3 object 'fwhm' with the following entries \cr
#' \itemize{
#'   \item fwhmInterpolator: A function that is used to find fwhm value for any given m/z.
#'   \item peaks: The m/z values of the peaks of a single spectrum that were used for fwhm calculation.
#'   \item fwhmVals: The corresponding fwhm values in Da.
#' }
#'
#' @export
#' @keywords internal
#'

fwhm <- function(fwhmInterpolator, peaks = NA_real_, fwhmVals = NA_real_) {

      # definition
      obj <- list(fwhmInterpolator = fwhmInterpolator, peaks = peaks, fwhmVals = fwhmVals)

      # integrity checks
      if(length(obj$peaks) != length(obj$fwhmVals)){
            stop("lengths of peaks and fwhmVals slots have to be equal. \n")
      }

      if(class(obj$fwhmInterpolator) != "function"){
            stop("provided fwhmInterpolator slot is not a function. \n")
      }


      class(obj) <- "fwhm"

      return(obj)

}


#' @export
print.fwhm <- function(obj) {

      cat("S3 object 'fwhm' holding fwhm data of a certain mass spectrum. \n ")
      cat("fwhm data is based on ", length(obj$peaks), "fwhm values,.\n")

}

#' @export
plot.fwhm <- function(obj, ...) {

   r = range(obj$peaks)
   qp = seq(r[1], r[2])
   plot(x = obj$peaks, y = obj$fwhmVals,
        main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm", ...)
   lines(x = qp, y = getFwhm(obj, qp), col = "green", lwd = 2)

}

#' Get fwhm for a given m/z value
#'
#' This is a method for S3 class `fwhm` to interpolate fwhm values.
#'
#' @param obj S3 object of type `fwhm`.
#' @param mzVals m/z value (Da) at which to estimate the fwhm.
#'
#' @return returns a numeric, estimated fwhm at the given m/z.
#'
#' @export
getFwhm <- function(obj, mzVals) {
   if(class(obj) != "fwhm")
      stop("obj must be an S3 object of type 'fwhm'\n")

   obj$fwhmInterpolator(mzVals)
}





#' sparseIntensityMatrix Class Constructor
#'
#' An S3 class to hold mass spectrometry imaging (MSI) intensity data in a sparse representation ('Matrix::dgCMatrix' class)
#' with few additional slots. It is used to speed-up computation. This constructor is for
#' internal use and normally should not be called by the user.
#'
#' @param spmat A sparse matrix of type 'Matrix::dgCMatrix' holding MSI intensity data.
#' @param mzAxis A numeric, the overall m/z axis for the experiment.
#' @param coordinates  A two-column data.frame, the coordinates of the spectra which are stored
#' in rows of 'spmat'.
#'
#' @return
#' An S3 object 'sparseIntensityMatrix' with the following entries \cr
#' \itemize{
#'   \item spmat: A sparse matrix of type 'Matrix::dgCMatrix' holding MSI intensity data.
#'   \item mzAxis: A numeric, the overall m/z axis for the experiment.
#'   \item coordinates: A two-column data.frame, the coordinates of the spectra which are stored in rows of 'spmat'.
#' }
#'
#' @export
#' @keywords internal
#'

sparseIntensityMatrix <- function(spmat, mzAxis, coordinates) {

# definition
obj <- list(spmat = spmat, mzAxis = mzAxis, coordinates = coordinates)

# integrity checks
if(class(obj$spmat) != "dgCMatrix"){
      stop("slot 'spmat' must be of type 'Matrix::dgCMatrix' \n")
}

if(length(obj$mzAxis) != obj$spmat@Dim[2]){
      stop("'mzAxis' slot length must be equal to the number of columns of the sparse data matrix 'spmat'.\n")
}

if(class(obj$coordinates) != "data.frame"){
      stop("slot 'coordinates' must be a data.frame. \n")
}

if(ncol(obj$coordinates) != 2){
      stop("slot 'coordinates' must be a data.frame of two columns. \n")
}

if(nrow(obj$coordinates) != obj$spmat@Dim[1]){
      stop("number of rows of @coordinates must be equal to the number of rows of the sparse matrix.\n")
}


class(obj) <- "sparseIntensityMatrix"

return(obj)

}

print.sparseIntensityMatrix <- function(obj) {

      cat("S3 object 'sparseIntensityMatrix' holding sparse intensity data of an MSI dataset. \n")
      cat("spmat: sparse matrix of type 'Matrix::dgCMatrix' with ", nrow(obj$spmat),
          "rows (spectra/pixels) and ", ncol(obj$spmat), "columns (m/z bins) .\n")
      cat("mzAxis: a numeric vector holding the m/z axis = number of columns of spmat.\n")
      cat("coordinates: a 2-column data.frame holding the spectra/pixel coordinates. \n")

}




#' analytePointPattern Class Constructor
#'
#' This is a class to extend 'ppp' class only enforcing specific point-marks. It also
#' adds a 'metaData' slot to declare which analytes are constructing a given 'analytePointPattern' object.
#' It represents the spatial point pattern of an analyte (or set of
#' analytes) with marks being a data.frame of two columns c('idx', 'intensity'). This constructor is for
#' internal use and normally should not be called by the user.
#'
#' @param spp a spatial point pattern of type 'ppp' with its `marks` being a dataframe with
#' two columns `idx`and `intensity`.
#' @param x,y,intensity,win optional, either these arguments or 'spp' have to be supplied.
#' @param metaData a named list having additional attributes describing the individual masses in `mzVals`.
#'
#' @details `mzVals` and `metaData` are combined into a data.frame which becomes a class slot `metaData`.
#' This contains the `mzVals` and possible other meta data with each row having a unique identifier `idx`
#' which then link `metaData` with `spp$marks`.
#'
#' @return
#' An S3 object of type 'ppp' and 'analytePointPattern' with the following entries \cr
#' \itemize{
#'   \item regular `ppp` entries: see `spatstat` documentation for `ppp` class.
#'   \item metaData: A data frame with meta data describing the points in `spp$marks`.
#' }
#'
#' @export
#' @keywords internal
#'
#'

analytePointPattern <- function(spp = NA, x = NA, y = NA, win = NA, intensity = NA, mzVals, metaData = list()) {

      #>>> to do: move mzVals to metaData > affected: searchAnalyte, superimposeAnalytes, ...

      # integrity checks
      if(all(identical(spp, NA), any(is.na(x), is.na(y), is.na(intensity), is.na(win)))){
         stop("either 'spp' or 'x,y,intensity,win' must be supplied. \n")
      }

      if(!is.list(metaData)){
         stop("metaData must be a named list. \n")
      }

      if(length(mzVals) < 1){
         stop("mzVals must be a numeric vector with at least one element. \n")
      }


      # definition
      idx <- as.numeric(Sys.time()) - prod(sample(seq(1,99), 3)) # primary key; id to connect to metaData table with marks (saves memory)

      if(identical(spp, NA)){


         spp <- ppp(x = x, y = y,
                              marks = data.frame(idx = rep(idx, length(intensity)), intensity = intensity),
                              window = win, checkdup = FALSE, drop = FALSE)


      } else {

         if(class(spp) != "ppp"){
            stop("slot 'spp' must be of type 'ppp' \n")
         }

         if(class(spp$marks) != "data.frame"){
            stop("marks of 'spp' must contain a data.frame of at least one column named 'intensity'.\n")
         }

         if(!("intensity" %in% colnames(spp$marks))){
            stop("marks of 'spp' must contain a data.frame of at least one column named 'intensity'.\n")
         }

         spp$marks <- data.frame(idx = rep(idx, length(spp$marks$intensity)), intensity = spp$marks$intensity)

      }

      # remove rejected points
      spp <- as.ppp(spp)

      # metaData object
      #mtdt <- data.frame(idx = idx, mzVals = mzVals, stringsAsFactors = FALSE)
      mtdt <- list(idx = idx, mzVals = mzVals)

      #if(length(metaData) > 0 & is.list(metaData)){
         #metaData <- metaData[!sapply(metaData, is.null)] # to remove residual null entries
         spp$metaData <- as.data.frame(c(mtdt, metaData),
                                       stringsAsFactors = FALSE)
      # } else {
      #    spp$metaData <- mtdt
      # }

      class(spp) <- c(class(spp), "analytePointPattern")

      return(spp)

}

#' Plot analytePointPattern
#'
#' A method to plot `anlaytePointPattern` objects. This calls `plot.ppp` with
#' pre-defined defaults. For more control over the plotting please use `plot.ppp`.
#'
#' @param obj S3 object of type `anlaytePointPattern` (and `spp`).
#' @param colourPal the colourmap to be used, see `?viridis::viridis_pal`.
#' @param uniformCol a character specifying a single colour. This will override `colourPal`.
#' @param transpFactor Transparency fraction. Numerical value or vector of values between 0 and 1,
#' giving the opaqueness of a colour. A fully opaque colour has `transpFactor=1`.
#' @param rescale logical, whether to scale the intensities of the output plot
#' to [0,1] interval.
#' @param pch     a positive integer, or a single character. See `?par`.
#' @param size    the size of the symbol: a positive number or zero. see`?symbolmap`.
#' @param main  character, title of the plot. If not given, the m/z value of `obj` is used.
#' @param leg.side position of legend relative to main plot.
#' @param leg.args list of additional arguments passed to control the legend. see
#' `?spatstat.geom::plot.ppp` for more details.
#' @param ... further arguments passed to `spatstat.geom::plot.ppp`.
#'
#' @return nothing, plots only.
#'
#' @export
#'
plotAnalyte <- function(obj, colourPal = "inferno", uniformCol = NULL, transpFactor = 0.7,
                        rescale = TRUE, pch = 19, size = 0.4, main = NULL,
                        leg.side = "right", leg.args = list(cex = 3, cex.axis = 1.25), ...){

  if(!("analytePointPattern" %in% class(obj))){
    stop("provided obj is not of type 'analytePointPattern'. \n")
  }

  if(is.null(main)){
        if(length(obj$metaData$mzVals) > 1){
              main <- paste0("Collective SPP of ", length(obj$metaData$mzVals)," m/z Values")
        } else{
                  mz <- ifelse(is.numeric(obj$metaData$mzVals),
                              as.character(round(obj$metaData$mzVals, 4)),
                              obj$metaData$mzVals) # to account for simulation cases

              main <- paste0("m/z ", mz)
        }

  }

   # msi uses reversed y axis
   args <- list(...)

   if("ylim" %in% names(args)){
    yrange <- rev(args[["ylim"]])
   } else {
    yrange <- rev(obj$window$yrange)
   }

  if(obj$n == 0){

    plot.ppp(obj, use.marks = TRUE, which.marks = "intensity",
             ylim = yrange,
             main = main)

  } else{

    if(is.null(uniformCol)){

             if(rescale){

                      obj <- .rescale(obj)
             }


      colfun <- colourmap(col = to.transparent((viridis::viridis_pal(option = colourPal)(100)), transpFactor),
                          range = range(obj$marks$intensity))

      plot.ppp(obj, use.marks = TRUE, which.marks = "intensity",
               ylim = yrange,
               main = main,
               symap = symbolmap(pch = pch,
                                 cols = colfun,
                                 size = size,
                                 range = range(obj$marks$intensity)),
               leg.side = leg.side,
               leg.args = leg.args,
               ...) # colors according to intensity

    } else{

             if(rescale){

                      obj <- .rescale(obj)
             }

      col <- to.transparent(uniformCol, transpFactor)

      plot.ppp(obj, use.marks = TRUE, which.marks = "intensity",
               ylim = yrange,
               main = main,
               symap = symbolmap(pch = pch,
                                 cols = col,
                                 size = size),
               leg.side = leg.side,
               leg.args = leg.args,
               ...)

    }


  }

}




#' Plot raster images of type `im`
#'
#' A method to plot `im` objects of the `spatstat` package. This calls `plot.im` with
#' pre-defined defaults. For more control over the plotting please use `plot.im`.
#'
#' @param obj S3 object of type `im`.
#' @param colourPal the colourmap to be used, see `?viridis::viridis_pal`.
#' @param uniformCol a character specifying a single colour. This will override `colourPal`.
#' @param rescale logical, whether to scale the intensities of the output plot
#' to [0,1] interval.
#' @param transpFactor Transparency fraction. Numerical value or vector of values between 0 and 1,
#' giving the opaqueness of a colour. A fully opaque colour has `transpFactor=1`.
#' @param irange a numeric of length 2, a custome intensity range for plotting. Incompatible with `rescale`;
#' if provided, `rescale` will be set to `FALSE`.
#' @param ribargs a list of additional arguments to control the display of the ribbon. see
#' `?spatstat.geom::plot.im` for details.
#' @param ribsep a numeric controlling the space between the ribbon and the image.
#' @param ... further arguments passed to `plot.im`.
#'
#' @return nothing, plots only.
#'
#' @export
#'
#'
plotImg <- function(obj, colourPal = "inferno", uniformCol = NULL,
                    rescale = TRUE, transpFactor = 0.7,
                    irange = NULL, ribargs = list(las = 2, cex.axis =1.25),
                    ribsep = 0.05, ...){

  if(class(obj) != "im"){
    stop("provided obj must be of type 'im'.\n")
  }

   # msi uses reversed y axis
   imageArgs <- list(...)

   if("ylim" %in% names(imageArgs)){
    yrange <- rev(imageArgs[["ylim"]])
   } else {
    yrange <- rev(obj$yrange)
   }


  if(is.null(uniformCol)){

    if(rescale & is.null(irange)){

      obj <- .rescale(obj)
    }

    if(is.null(irange)){
      irange <- range(obj$v, na.rm = T)
    }



    plot.im(x = obj,
            col = colourmap(viridis::viridis_pal(option = colourPal)(100),
                            range = irange),
            ylim = yrange,
            box = FALSE,
            ribargs = ribargs,
            ribsep = ribsep,
            ...)


  } else{



    col <- to.transparent(uniformCol, transpFactor)


    plot.im(x = obj,
            col = col,
            ylim = yrange,
            box = FALSE,
            ribargs = ribargs,
            ribsep = ribsep,
            ...)

  }



}


#' Convert `anlaytePointPattern` to `im`
#'
#' A method to convert `anlaytePointPattern` objects to `im` objects.
#'
#' @param obj S3 object of type `anlaytePointPattern` (and `spp`).
#' @param rescale logical, whether to scale the intensities of the output `im` object
#' to [0,1] interval.
#' @param zero.rm for internal use only.
#' @param ... additional arguments passed to `spatstat.geom::pixellate`.
#'
#' @return an object of type `im` (see `spatstat` documentation).
#'
#' @export
#' @keywords internal
#'

spp2im <- function(obj, rescale = FALSE, zero.rm = FALSE, ...){

         if(!("analytePointPattern" %in% class(obj))){
                  stop("provided obj is not of type 'analytePointPattern'. \n")
         }

         win <- obj$window

         im <- pixellate(obj,
                         weights = obj$marks$intensity,
                         W = as.mask(win,dimyx=c(diff(win$yrange) + 1, diff(win$xrange) + 1)),
                         ...)

         im <- as.im(im) # fix coordinates inconsistancies

         if(zero.rm){
                  im[im == 0] <- NA
         }


         if(rescale){

                  im <- .rescale(im)
         }

         return(im)
}

.rescale <- function(obj){

         if(identical(class(obj), "im")){

                  im <- (obj - min(obj)) / (max(obj) - min(obj))

                  return(im)
         }

         if("analytePointPattern" %in% class(obj)){

                  intens <- obj$marks$intensity

                  obj$marks$intensity <- (intens - min(intens)) / (max(intens) - min(intens))

                  return(obj)

         }
}


#' molProbMap Class Constructor
#'
#' This is an S3 class to represent the calculated molecular probability maps. This
#' class is only used with objects returned by `porbMap` function and should not be
#' called explicitly by the user.
#'
#' @param bw   A numeric, the calculated Gaussian bandwidth
#' @param sppMoi A 'analytePointPattern' object representing a given analyte (metabolite-of-interest).
#' @param csrMoi A 'analytePointPattern' object representing the corresponding CSR model.
#' @param rhoMoi An 'im' object storing the density image of 'sppMoi' after applying KDE.
#' @param rhoCsr An 'im' object storing the density image of 'csrMoi' after applying KDE.
#' @param hotspotpp A 'analytePointPattern' object holding only the points that are located within the hotspot.
#' @param hotspotIm An 'im' object storing the corresponding pixellated image of the 'hotspotpp' slot.
#' @param hotspotMask An 'owin' object storing the corresponding window of the 'hotspotpp' slot.
#' @param coldspotpp A 'analytePointPattern' object holding only the points that are located within the coldspot.
#' @param coldspotIm An 'im' object storing the corresponding pixellated image of the 'coldspotpp' slot.
#' @param coldspotMask An 'owin' object storing the corresponding window of the 'coldspotpp' slot.
#'
#' @return
#' An S3 object of type 'molProbMap'.
#'
#' @export
#' @keywords internal
#'
molProbMap <- function(bw, sppMoi, csrMoi, rhoMoi, rhoCsr, hotspotpp, hotspotIm, hotspotMask, coldspotpp, coldspotIm, coldspotMask) {

      # definition
      obj <- list(bw = bw, sppMoi = sppMoi, csrMoi = csrMoi, rhoMoi = rhoMoi, rhoCsr = rhoCsr,
                  hotspotpp = hotspotpp, hotspotIm = hotspotIm, hotspotMask = hotspotMask,
                  coldspotpp = coldspotpp, coldspotIm = coldspotIm, coldspotMask = coldspotMask)

      # integrity checks
      if(!("analytePointPattern" %in% class(obj$sppMoi))){
            stop("slot 'sppMoi' must be of type 'ppp' and 'analytePointPattern' \n")
      }

      if(!("analytePointPattern" %in% class(obj$csrMoi))){
         stop("slot 'csrMoi' must be of type 'analytePointPattern' \n")
      }

      if(class(obj$rhoMoi) != "im"){
         stop("slot 'rhoMoi' must be of type 'im' from 'spaststat' \n")
      }

      if(class(obj$rhoCsr) != "im"){
         stop("slot 'rhoCsr' must be of type 'im' from 'spaststat' \n")
      }

      if(!("analytePointPattern" %in% class(obj$hotspotpp))){
         stop("slot 'hotspotpp' must be of type 'analytePointPattern' \n")
      }

      if(class(obj$hotspotIm) != "im"){
         stop("slot 'hotspotIm' must be of type 'im' from 'spaststat' \n")
      }

      if(class(obj$hotspotMask) != "owin"){
         stop("slot 'hotspotMask' must be of type 'owin' from 'spaststat' \n")
      }

      if(!("analytePointPattern" %in% class(obj$coldspotpp))){
         stop("slot 'coldspotpp' must be of type 'analytePointPattern' \n")
      }

      if(class(obj$coldspotIm) != "im"){
         stop("slot 'coldspotIm' must be of type 'im' from 'spaststat' \n")
      }

      if(class(obj$coldspotMask) != "owin"){
         stop("slot 'coldspotMask' must be of type 'owin' from 'spaststat' \n")
      }

      class(obj) <- "molProbMap"

      return(obj)

}

#' @export
print.molProbMap <- function(obj) {

      cat("S3 object 'molProbMap' holding molecular probability maps of a given analyte in an MSI dataset. \n ")

}


#' Plot molecular Probability Maps
#'
#' A method to plot MPMs with different combinations.
#'
#' @param obj S3 object of type `molProbMap`.
#' @param what What to plot, c("detailed", "analytePointPattern", "csrPointPattern", "analyteDensityImage",
#' "csrDensityImage", "kdeIntensitiesDistr", "MPM"). The "detailed" option plots all options.
#' @param sppArgs a named list of arguments to be passed to `plotAnalyte` or plotting
#' "analytePointPattern" and "csrPointPattern". There are sensible defaults.
#' @param imageArgs a named list of arguments to be passed to `spatstat.geom::plot.im` for plotting
#' "analyteDensityImage", "csrDensityImage" and "MPM". There are sensible defaults.
#' @param fArgs a named list of arguments to be passed to `base::plot` for plotting
#' "kdeIntensitiesDistr". There are sensible defaults.
#' @param rescale logical, whether to scale the intensities of the output plot
#' to [0,1] interval.
#' @param signifArea a character indicating which significance area to plot i.e. c("both", "hotspot", "coldspot").
#' @param ionImage an optional rastered image of type `im` of the corresponding "regular"
#' ion image, used for comparison. Could be generated via `moleculaR::searchAnalyte(..., wMethod = "sum")`
#' and subsequently using `pixellate`.
#' @param figLegend logical, whether to show the hotspot and coldspot legends.
#' @param lwd.signifArea numeric, the width of the hotspot and coldspot contours.
#' @param col.hotspot a character specifying the color of the hotspot contour.
#' @param col.coldspot a character specifying the color of the coldspot contour.
#'
#' @return nothing, plots only.
#'
#' @method plot molProbMap
#'
#' @export
#'
plot.molProbMap <- function(obj, what = "detailed",
                            sppArgs = list(),
                            imageArgs = list(),
                            fArgs = list(),
                            rescale = TRUE,
                            signifArea = "both",
                            ionImage = NA,
                            figLegend = TRUE,
                            lwd.signifArea = 5,
                            col.hotspot = "red",
                            col.coldspot = "blue"){

         # checks
         if(!(signifArea %in% c("hotspot", "coldspot", "both"))){
                  stop("'signifArea' must be one of c('hotspot', 'coldspot', 'both').\n")
         }



   switch (what,
      "analytePointPattern" = {

               if(rescale){
                        obj$sppMoi <- .rescale(obj$sppMoi)
               }

               # workout the plotting arguments
               # sppDefaults    <- list(colourPal = "inferno", uniformCol = NULL, transpFactor = 0.7,
               #                        pch = 19, size = 0.4, main = NULL)
               # sppArgs        <- .mergeArgs(sppArgs, sppDefaults)

               # call the plotting function
               do.call(plotAnalyte, c(list(obj = obj$sppMoi), sppArgs))


      },
      "csrPointPattern" = {

               if(rescale){
                        obj$csrMoi <- .rescale(obj$csrMoi)
               }

               # workout the plotting arguments
               # sppDefaults    <- list(colourPal = "inferno", uniformCol = NULL, transpFactor = 0.7,
               #                        pch = 19, size = 0.4, main = NULL)
               # sppArgs        <- .mergeArgs(sppArgs, sppDefaults)



               # call the plotting function
               do.call(plotAnalyte, c(list(obj = obj$csrMoi), sppArgs))


      },
      "analyteDensityImage" = {

               if(rescale){
                        obj$rhoMoi <- .rescale(obj$rhoMoi)
               }

               #workout the plotting arguments
               imageDefaults  <- list(main = expression(paste(rho["MOI"], "(x,y)")))

               imageArgs      <- .mergeArgs(imageArgs, imageDefaults)


               # call the plotting function
               do.call(plotImg, c(list(obj = obj$rhoMoi), imageArgs))


      },
      "csrDensityImage" = {

               if(rescale){
                        obj$rhoCsr <- .rescale(obj$rhoCsr)
               }

               # workout the plotting arguments
               imageDefaults  <- list(main = expression(paste(rho["CSR"], "(x,y)")))

               imageArgs      <- .mergeArgs(imageArgs, imageDefaults)




               # call the plotting function
               do.call(plotImg, c(list(obj = obj$rhoCsr), imageArgs))


      },
      "kdeIntensitiesDistr" = {

               if(rescale){
                        obj$rhoMoi <- .rescale(obj$rhoMoi)
                        obj$rhoCsr <- .rescale(obj$rhoCsr)
               }

               fmoi <- density(c(obj$rhoMoi$v), na.rm = TRUE)
               fcsr <- density(c(obj$rhoCsr$v), na.rm = TRUE)
               muCsr <- mean(obj$rhoCsr$v, na.rm = TRUE)
               sigmaCsr <- sd(obj$rhoCsr$v, na.rm = TRUE)
               lowerCutoff <- qnorm(0.05, muCsr, sigmaCsr)
               upperCutoff <- qnorm(0.05, muCsr, sigmaCsr, lower.tail = FALSE)

               ymax <- max(max(fmoi$y), max(fcsr$y))
               ymin <- min(min(fmoi$y), min(fcsr$y))
               xmax <- max(max(fmoi$x), max(fcsr$x))
               xmin <- min(min(fmoi$x), min(fcsr$x))

               # Create two panels side by side
               #layout(t(1:1), widths=c(5,1))

               # Set margins and turn all axis labels horizontally (with `las=1`)
               #par(mar=rep(3, 4), oma=rep(4, 4), las=1)

               # workout the plotting arguments
               fDefaults  <- list(main = bquote(f[CSR]*(k) ~ and ~ f[MOI]*(k) ~ at ~ bw==.(obj$bw)),
                                  col = "black",
                                  ylim = c(ymin, ymax),
                                  xlim = c(xmin, xmax),
                                  xlab = "Intensities",
                                  ylab = "Density",
                                  bty = "n",
                                  lwd = 3)

               fArgs      <- .mergeArgs(fArgs, fDefaults)

               # call the plotting function
               do.call(plot, c(list(fcsr), fArgs))


               lines(fmoi, col =  "forestgreen", lwd = fArgs$lwd)

               polygon(x = c(upperCutoff, xmax, xmax, upperCutoff),
                       y = c(ymin, ymin, ymax, ymax), col = to.transparent("coral2", 0.5),
                       border = NA)

               polygon(x = c(xmin, lowerCutoff, lowerCutoff, xmin),
                       y = c(ymin, ymin, ymax, ymax), col = to.transparent("cornflowerblue", 0.5),
                       border = NA)

               polygon(x = c(fcsr$x[fcsr$x >= upperCutoff],  upperCutoff),
                       y = c(fcsr$y[fcsr$x >= upperCutoff],  0),
                       col = to.transparent("red", 0.5),
                       border = NA)

               polygon(x = c(fcsr$x[fcsr$x <= lowerCutoff], lowerCutoff),
                       y = c(fcsr$y[fcsr$x <= lowerCutoff],  0),
                       col = to.transparent("blue", 0.5),
                       border = NA)

               if(figLegend){
                        legend("topright", legend = c(expression(paste(italic(f["CSR"]), italic("(k)"))),
                                                      expression(paste(italic(f["MOI"]), italic("(k)")))),
                               lty = c("solid"), lwd = 2,
                               col = c("black", "forestgreen"),
                               bty = "n", horiz = FALSE)
               }


      },
      "MPM" = {

         #spwin <- obj$sppMoi$window

         if(length(obj$sppMoi$metaData$mzVals) > 1 | obj$sppMoi$metaData$mzVals[1] == "expr"){
            # this indicates CPPM -> use the density image

            imgMpm <- obj$rhoMoi

         } else { # this is a single-analyte MPM -> use rasterized image

            # raster image of the spp
            # imgMpm  <- pixellate(obj$sppMoi,
            #                      weights = obj$sppMoi$marks$intensity,
            #                      W = as.mask(spwin,dimyx=c(diff(spwin$yrange) + 1, diff(spwin$xrange) + 1)),
            #                      padzero = FALSE, savemap = FALSE)

            imgMpm <- spp2im(obj$sppMoi)

         }

         if(rescale){
                  imgMpm <- .rescale(imgMpm)
         }

         # workout the plotting arguments
         imageDefaults  <- list(main = NULL)

         imageArgs      <- .mergeArgs(imageArgs, imageDefaults)

         # fix title
         if(is.null(imageArgs$main)){
                  if(length(obj$sppMoi$metaData$mzVals) > 1){
                           imageArgs$main <- paste0("CPPM of ", length(obj$sppMoi$metaData$mzVals)," m/z Values")
                  } else{
                           mz <- ifelse(is.numeric(obj$sppMoi$metaData$mzVals),
                                        as.character(round(obj$sppMoi$metaData$mzVals, 4)),
                                        obj$sppMoi$metaData$mzVals) # to account for simulation cases

                           imageArgs$main <- paste0("m/z ", mz, " MPM")
                  }
         }



         # call the plotting function
         do.call(plotImg, c(list(obj = imgMpm), imageArgs))


         if(signifArea == "both" | signifArea == "hotspot"){
            plot.owin(obj$hotspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = lwd.signifArea,  add = TRUE)
            plot.owin(obj$hotspotpp$window, col = rgb(1,1,1,0.0), border = col.hotspot, lwd = lwd.signifArea/2, lty = "dashed",add = TRUE)
         }

         if(signifArea == "both" | signifArea == "coldspot"){
            plot.owin(obj$coldspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = lwd.signifArea,  add = TRUE)
            plot.owin(obj$coldspotpp$window, col = rgb(1,1,1,0.0), border = col.coldspot, lwd = lwd.signifArea/2, lty = "dashed",add = TRUE)
         }

         if(figLegend){
                  legend("bottom", legend = c("Analyte Hotspot", "Analyte Coldspot"), lty = c("dashed"),
                           col = c("red", "blue"), bty = "n", horiz = TRUE, inset = c(-0.5,0))

         }


      },
      "detailed" = {

         par(mfrow = c(3, 2))

         plot(obj = obj, what = "csrPointPattern",
              sppArgs = sppArgs,
              imageArgs = imageArgs,
              fArgs = fArgs,
              rescale = rescale,
              signifArea = signifArea,
              ionImage = ionImage,
              figLegend = figLegend,
              lwd.signifArea = lwd.signifArea,
              col.hotspot = col.hotspot,
              col.coldspot = col.coldspot)

         plot(obj = obj, what = "analytePointPattern",
              sppArgs = sppArgs,
              imageArgs = imageArgs,
              fArgs = fArgs,
              rescale = rescale,
              signifArea = signifArea,
              ionImage = ionImage,
              figLegend = figLegend,
              lwd.signifArea = lwd.signifArea,
              col.hotspot = col.hotspot,
              col.coldspot = col.coldspot)

         plot(obj = obj, what = "csrDensityImage",
              sppArgs = sppArgs,
              imageArgs = imageArgs,
              fArgs = fArgs,
              rescale = rescale,
              signifArea = signifArea,
              ionImage = ionImage,
              figLegend = figLegend,
              lwd.signifArea = lwd.signifArea,
              col.hotspot = col.hotspot,
              col.coldspot = col.coldspot)

         plot(obj = obj, what = "analyteDensityImage",
              sppArgs = sppArgs,
              imageArgs = imageArgs,
              fArgs = fArgs,
              rescale = rescale,
              signifArea = signifArea,
              ionImage = ionImage,
              figLegend = figLegend,
              lwd.signifArea = lwd.signifArea,
              col.hotspot = col.hotspot,
              col.coldspot = col.coldspot)

         if(identical(ionImage, NA)){
                  plot(obj = obj, what = "kdeIntensitiesDistr",
                       sppArgs = sppArgs,
                       imageArgs = imageArgs,
                       fArgs = fArgs,
                       rescale = rescale,
                       signifArea = signifArea,
                       ionImage = ionImage,
                       figLegend = figLegend,
                       lwd.signifArea = lwd.signifArea,
                       col.hotspot = col.hotspot,
                       col.coldspot = col.coldspot)
         }

         plot(obj = obj, what = "MPM",
              sppArgs = sppArgs,
              imageArgs = imageArgs,
              fArgs = fArgs,
              rescale = rescale,
              signifArea = signifArea,
              ionImage = ionImage,
              figLegend = figLegend,
              lwd.signifArea = lwd.signifArea,
              col.hotspot = col.hotspot,
              col.coldspot = col.coldspot)

         if(!(identical(ionImage, NA))){

                  if(rescale){
                           ionImage <- .rescale(ionImage)
                  }

                  # workout the plotting arguments
                  imageDefaults  <- list(main = "Corresponding Ion Image")

                  imageArgs      <- .mergeArgs(imageArgs, imageDefaults)



                  # call the plotting function
                  do.call(plotImg, c(list(obj = ionImage), imageArgs))

         }


      },
      stop("argument 'what' must be one of c('detailed', 'analytePointPattern', 'csrPointPattern', ",
      "'analyteDensityImage','csrDensityImage', 'kdeIntensitiesDistr' or 'MPM').\n")
   )

}

# merge default arguments with provided ones
# providedArgs and defualt args are named lists
.mergeArgs <- function(providedArgs, defaultArgs){

         if(length(providedArgs) == 0){
                  providedArgs <- defaultArgs
         } else {
                  providedArgs <- c(providedArgs, defaultArgs[!(names(defaultArgs) %in% names(providedArgs))])
         }

         return(providedArgs)
}

# .assignTitle <- function(mzVals, expr = NULL) {
#
#          if(!is.null(expr)){
#                   return(expr)
#          }
#
#
#          if(length(mzVals) > 1){
#                   title <- paste0("Collective SPP of ", length(mzVals)," m/z Values")
#          } else{
#                   mz <- ifelse(is.numeric(mzVals),
#                                as.character(round(mzVals, 4)),
#                                mzVals) # to account for simulation cases
#
#                   title <- paste0("m/z ", mz, " SPP")
#          }
#
#          return(title)
# }

#' lipidSearchList Class Constructor
#'
#' A class to re-arrange lipid listings in the SwissLipid database to facilitate computational
#' serach. It splits the full database into lipid class-oriented lists. This function is for
#' internal use and normally should not be called by the user.
#'
#' @param lipidList a data.frame which lists entries of lipids.
#' @param hitsList a data.frame holding the lipids that were detected.
#' @param allClasses a character vector summarizing all classes taking part in the lipid search.
#'
#' @return
#' An S3 object of type 'lipidSearchList'.
#'
#' @export
#' @keywords internal
#'

lipidSearchList <- function(lipidList, hitsList, allClasses) {


   # definition
   l <- list(lipidList = lipidList, hitsList = hitsList, allClasses = allClasses)
   class(l) <- "lipidSearchList"

   return(l)

}
