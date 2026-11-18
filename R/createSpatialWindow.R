#' Create Spatial Window of MSI Data
#'
#' This function creates `spatstat.geom::owin` polygonal spatial window which represents
#' the boundaries of the tissue under study.
#'
#' @param pixelCoords: 	a data frame with x and y columns containing the coordinates
#' of all pixels within a given MSI data. A `matrix` will be coerced into a
#' data frame.
#' @param isPolygon: a logical, whether `pixelCoords` represents the coordinates
#' of a polygon or not.
#' @param clean: logical, whether to fill holes created by missing pixels.
#' @param plot:  whether to plot.
#' @param pixelSize Optional positive pixel width (one number for square pixels,
#' or x and y widths). For pixel coordinates, fixes the grid spacing instead of
#' inferring it from observed positions. Supports single rows or columns and
#' gaps between observed pixels. Ignored when \code{isPolygon = TRUE}.
#' @return
#' A spatial window object of type `spatstat.geom::owin`.
#'
#' @export
#'
createSpatialWindow <- function(pixelCoords,
                                isPolygon = FALSE,
                                clean = FALSE,
                                plot = TRUE, pixelSize = NULL) {


      #// checks
      if(!is.data.frame(pixelCoords)){
            if(is.matrix(pixelCoords)){
                  pixelCoords <- as.data.frame(pixelCoords)
            } else {
                  stop("pixelCoords must either be a matrix or data.frame ",
                       "with x and y columns containing pixel coordinates. \n")
            }

      }

      if(ncol(pixelCoords) != 2){
            stop("pixelCoords must be a matrix or data.frame of two columns named 'x' and 'y'. \n")
      }

      # create the spwin from the raster
      if(isPolygon){

            if(!spatstat.utils::verify.xypolygon(pixelCoords)){
                  stop("PixelCoords is not a polygon. \n")
            }

            if(spatstat.utils::Area.xypolygon(pixelCoords) < 0){
                  warning("The plogon provided as pixelCoords has a negative area i.e. ",
                          "traversed clockwise (which for spatstat represents a hole). ",
                          "The order of the points was reversed. \n")

                  pixelCoords <- spatstat.utils::reverse.xypolygon(pixelCoords)

            }

            spwin       <- spatstat.geom::owin(poly = pixelCoords)

      } else {


            if(is.null(pixelSize)){
                  spwin <- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = pixelCoords))
            } else {
                  if(!is.numeric(pixelSize) || !(length(pixelSize) %in% c(1L, 2L)) ||
                     any(!is.finite(pixelSize)) || any(pixelSize <= 0)){
                        stop("pixelSize must contain one or two positive finite numbers.")
                  }
                  pixelSize <- rep(pixelSize, length.out = 2L)
                  coords <- as.matrix(pixelCoords)
                  if(!is.numeric(coords) || nrow(coords) == 0L || any(!is.finite(coords))){
                        stop("pixelCoords must contain finite numeric pixel coordinates.")
                  }
                  origin <- apply(coords, 2, min)
                  offsets <- sweep(sweep(coords, 2, origin), 2, pixelSize, "/")
                  if(any(abs(offsets - round(offsets)) > 1e-7)){
                        stop("pixelCoords do not lie on the grid specified by pixelSize.")
                  }
                  indices <- round(offsets) + 1L
                  extent <- apply(indices, 2, max)
                  # An empty border lets spatstat polygonize single rows/columns.
                  mask <- matrix(FALSE, nrow = extent[2] + 2L, ncol = extent[1] + 2L)
                  mask[indices[, 2:1, drop = FALSE] + 1L] <- TRUE
                  spwin <- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = mask,
                        xrange = origin[1] + c(-1.5, extent[1] + 0.5) * pixelSize[1],
                        yrange = origin[2] + c(-1.5, extent[2] + 0.5) * pixelSize[2]))
            }

      }




      # remove holes
      if(clean){

            owinInfo <- summary(spwin)

            tokp <- vector("logical", owinInfo$npoly)

            for(i in 1 : length(tokp)){

                  tokp[i] <- ifelse(spatstat.utils::Area.xypolygon(spwin$bdry[[i]]) > 0, TRUE, FALSE)

            }

            spwin <- spatstat.geom::owin(poly = spwin$bdry[tokp])

      }




      if(plot) {


            spatstat.geom::plot.owin(spwin, ylim = rev(spwin$yrange), main = "Spatial Window")


      }



      return(spwin)




}



