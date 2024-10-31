
#' @include manualSpatstatImport.R
NULL


#' Normalized cross correlation
#'
#' This function Calculates the normalized cross correlation between two
#' images of type `im` (see `spatstat` documentation).
#'
#' @param im0:    a reference image object of type `im`.
#' @param im1:    a second image object of type `im`.
#' @return
#' A numeric, the calculated normalized cross correlation
#'
#' @export
#' @keywords internal
#'
ncc <- function(im0, im1){

      stopifnot(class(im0) == "im", class(im1) == "im")

      # have to have the same scale
      im0 <- .rescale(im0)
      im1 <- .rescale(im1)

      mu0 <- mean(im0)
      mu1 <- mean(im1)

      num <- sum((im0 - mu0) * (im1 - mu1))
      denom <- sqrt(sum((im0-mu0)^2) * sum((im1-mu1)^2))

      return(num/denom)

}


#' Mean Squared Error
#'
#' This function Calculates the mean squared error
#'  between two
#' images of type `im` (see `spatstat` documentation).
#'
#' @param im0:    a reference image object of type `im`.
#' @param im1:    a second image object of type `im`.
#' @return
#' A numeric, the calculated mean squared error.
#'
#' @export
#' @keywords internal
#'
mse <- function(im0, im1){

      stopifnot(class(im0) == "im", class(im1) == "im")

      # have to have the same scale
      im0 <- .rescale(im0)
      im1 <- .rescale(im1)

      val <- sum((im0 - im1)^2)/prod(dim(im0))

      return(val)



}


#' Dice Similarity Coefficient of two window objects
#'
#' This function Calculates the Dice similarity coefficient of two
#' window objects of type `owin` (see `?spatstat.geom::owin`).
#'
#' @param win0:    a reference image object of type `owin`.
#' @param win1:    a second image object of type `owin`.
#' @param winBackground: an optional bounding window of type `owin` which acts as
#' the background window for both `win0` and `win1`. This is only used for
#' visualization.
#' @param plot:   whether to plot the result.
#' @return
#' A numeric, the calculated Dice similarity coefficient.
#'
#' @export
#' @keywords internal
#'
dsc <- function(win0, win1, winBackground = NULL, plot = FALSE){

      stopifnot(class(win0) == "owin", class(win1) == "owin")

      s   <- intersect.owin(win0, win1)

      if(is.null(s)) {

            warning("No intersection of win0 and win1, zero DSC is returned.\n")

            return(0)

      }

      dscVal   = (2 * area.owin(s))/(area.owin(win0) + area.owin(win1))

      if(plot){
            if(!is.null(winBackground)){
                  plot.owin(winBackground,
                            ylim = rev(winBackground$yrange),
                            main = paste0("Dice Similarity Coefficient = ", round(dscVal, 3)))
                  plot.owin(win0, col = rgb(0,1,0,1), add = TRUE)
                  plot.owin(win1, col = rgb(1,0,0,1), add = TRUE)
                  plot.owin(s, col = rgb(1,1,0,1), add = TRUE)
            } else {
                  plot.owin(win0,
                            ylim = rev(win0$yrange),
                            col = rgb(0,1,0,1),
                            main = paste0("Dice Similarity Coefficient = ", round(dscVal, 3)))
                  plot.owin(win1, col = rgb(1,0,0,1), add = TRUE)
                  plot.owin(s, col = rgb(1,1,0,1), add = TRUE)
            }

            legend("topleft", bty = "n", horiz = FALSE,
                   legend = c("win0", "win1", "intersection"),
                   col = c( rgb(0,1,0,1), rgb(1,0,0,1), rgb(1,1,0,1)),
                   pch = 15, cex = 1)



      }

      return(dscVal)



}



#' Convert `anlaytePointPattern` to `im`
#'
#' A method to convert `anlaytePointPattern` objects to `im` objects.
#'
#' @param obj S3 object of type `anlaytePointPattern` (and `spp`).
#' @param weighted a logical, whether to scale points by their intensities.
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

spp2im <- function(obj, weighted = TRUE, rescale = FALSE, zero.rm = FALSE, ...){

      if(!("analytePointPattern" %in% class(obj))){
            stop("provided obj is not of type 'analytePointPattern'. \n")
      }

      win <- obj$window

      if(weighted){
            w <- obj$marks$intensity
      } else{
            w <- NULL
      }

      im <- pixellate(obj,
                      weights = w,
                      #W = as.mask(win,dimyx=c(diff(win$yrange) + 1, diff(win$xrange) + 1)),
                      W = as.mask(win,dimyx=c(diff(win$yrange), diff(win$xrange))),
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

#' Convert `im` to `anlaytePointPattern`
#'
#' A method to convert `im` objects to `anlaytePointPattern` objects.
#'
#' @param obj S3 object of type `im`.
#' @param win S3 object of type `owin` representing the window or the resulting
#' `anlaytePointPattern` object.
#' @param rescale logical, whether to scale the intensities of the output `im` object
#' to [0,1] interval.
#' @param zero.rm for internal use only.
#' @param mzVal optional, the m/z value of the analyte represnted in `obj`.
#'
#' @return an object of type `anlaytePointPattern`.
#'
#' @export
#' @keywords internal
#'

im2spp <- function(obj, win , rescale = FALSE, zero.rm = TRUE, mzVal = "mz"){

      if(!(class(obj) == "im")){
            stop("provided obj is not of type 'im'. \n")
      }


      # remove zeroes which might be a result of ceorcing from ppp to im
      if(zero.rm){
            obj[obj == 0] <- NA
      }

      obj <- as.im(obj) # fix coordinates inconsistancies



      if(rescale){

            obj <- .rescale(obj)
      }

      # convert back to dataframe
      imdf <- as.data.frame.im(x = obj)

      if(any(is.infinite(imdf$value))){
            imdf <- imdf[which(is.finite(imdf$value)), ]
      }


      sppObj <- analytePointPattern(x = imdf$x, y = imdf$y,
                                    intensity = imdf$value,
                                    win = win, mzVals = "expr")


      return(sppObj)
}



.rescale <- function(obj, from = 0, to = 10){

      if(identical(class(obj), "im")){

            # offset the lowest value (zero) because it is reserved to pixels with empty points in spatstat
            minVal <- min(obj) - 0.0001

            im <- ((obj - minVal) / (max(obj) - minVal)) * (to - from) + from

            return(im)
      }

      if("analytePointPattern" %in% class(obj)){

            intens <- obj$marks$intensity

            # offset the lowest value (zero) because it is reserved to pixels with empty points in spatstat
            minVal <- min(intens) - 0.0001

            obj$marks$intensity <- ((intens - minVal) / (max(intens) - minVal)) * (to - from) + from

            return(obj)

      }
}


