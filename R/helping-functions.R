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
