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
#' @return
#' A spatial window object of type `spatstat.geom::owin`.
#'
#' @export
#'
createSpatialWindow <- function(pixelCoords,
                                isPolygon = FALSE,
                                clean = FALSE,
                                plot = TRUE) {


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


            spwin       <- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = pixelCoords))

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



