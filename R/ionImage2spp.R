#' Convert an ion image matrix or array to an analyte point pattern
#'
#' Convert a plain numeric intensity matrix from any image provider to the
#' spatial representation accepted by [probMap()]. Three-dimensional arrays,
#' including output from `MALDIquant::msiSlices()`, are also supported.
#'
#' @param x A numeric matrix or a three-dimensional numeric array. Missing
#'   pixels outside the tissue must be `NA` (or `NaN`).
#' @param mzVal A single m/z value or descriptive label. If `NULL`, use the
#'   selected entry of the array's `center` attribute (as supplied by
#'   `MALDIquant::msiSlices()`), or the label `"image"` if absent.
#' @param slice Index along the third dimension. Required for arrays with
#'   multiple layers; a singleton third dimension is selected automatically.
#'   Leave `NULL` for a plain matrix.
#'   The third dimension represents mass layers, not spatial z coordinates.
#' @param orientation Dimension convention: `"xy"` (default) maps the first
#'   dimension to x and the second to y, matching `MALDIquant::msiSlices()`.
#'   Use `"yx"` for matrices whose rows represent y and columns represent x.
#' @param win Optional tissue window of class `spatstat.geom::owin`. It must
#'   contain all non-missing pixel centres. Converted to a polygonal window.
#'   If `NULL`, build a window from all non-missing pixels with
#'   [createSpatialWindow()], preserving holes and disconnected regions.
#' @param zero.rm Logical; omit zero-intensity points (default `TRUE`), matching
#'   [searchAnalyte()]. Zero-valued pixels still contribute to the tissue window.
#'
#' @details
#' For a plain matrix, no mass metadata or spectrum input is required.
#' `mzVal` is optional and defaults to the descriptive label `"image"` when
#' there is no `center` attribute. Use `orientation = "yx"` for a conventional
#' matrix with rows representing y and columns representing x; the default
#' `"xy"` instead treats rows as x and columns as y.
#'
#' Coordinates are one-based array indices with unit pixel spacing. Dimension
#' names are not interpreted as coordinates. No transpose, vertical flip,
#' normalization, smoothing, or FWHM-dependent mass weighting is applied.
#' The resulting MPM uses the input image's existing intensities and mass-window
#' aggregation, which may differ from [searchAnalyte()]'s Gaussian weighting.
#'
#' Only the selected layer is converted. All non-missing intensities in that
#' layer must be finite and non-negative, with at least one positive value.
#' An all-missing or entirely zero image raises an error.
#' A window inferred from an image is appropriate only when its missing pixels
#' describe the tissue boundary. Supply a common `win` when tissue boundaries
#' are known separately or several images should share the same analysis domain.
#'
#' @return An `analytePointPattern` (also a `spatstat` `ppp`) containing original
#'   intensities, pixel coordinates, a polygonal window, and m/z metadata.
#' @seealso [im2spp()] for `spatstat` images, [probMap()], [createSpatialWindow()].
#' @examples
#' # A plain matrix: rows = y, columns = x.
#' image <- matrix(NA_real_, nrow = 20, ncol = 15)
#' image[3:18, 3:13] <- 1
#' image[7:13, 6:10] <- 10
#' spp <- ionImage2spp(image, orientation = "yx")
#' mpm <- probMap(spp, bw = 1, seed = 1, verbose = FALSE)
#'
#' # Optionally attach an m/z value or descriptive label.
#' spp <- ionImage2spp(image, orientation = "yx", mzVal = 788.5447)
#'
#' # Arrays are also supported; choose a layer when there is more than one.
#' layers <- array(rep(image, 2), dim = c(20, 15, 2))
#' attr(layers, "center") <- c(788.5447, 810.6)
#' spp <- ionImage2spp(layers, orientation = "yx", slice = 2)
#' # MALDIquant::msiSlices() uses [x, y, mass]: keep orientation = "xy" for it.
#' @export
#' @md
ionImage2spp <- function(x, mzVal = NULL, slice = NULL,
                          orientation = c("xy", "yx"), win = NULL,
                          zero.rm = TRUE) {
      dimensions <- dim(x)
      if(!is.numeric(x) || !(length(dimensions) %in% c(2L, 3L)) ||
         any(dimensions == 0L)){
            stop("x must be a non-empty numeric matrix or three-dimensional array.")
      }
      orientation <- match.arg(orientation)
      if(!is.logical(zero.rm) || length(zero.rm) != 1L || is.na(zero.rm)){
            stop("zero.rm must be TRUE or FALSE.")
      }
      centres <- attr(x, "center", exact = TRUE)
      layer <- 1L
      if(length(dimensions) == 3L){
            if(is.null(slice)){
                  if(dimensions[3] != 1L){
                        stop("x has multiple layers; select one explicitly with slice.")
                  }
                  slice <- 1L
            }
            if(!is.numeric(slice) || length(slice) != 1L || !is.finite(slice) ||
               slice < 1 || slice > dimensions[3] || slice != floor(slice)){
                  stop("slice must be a single integer index within the third dimension.")
            }
            layer <- slice
            # Preserve both spatial dimensions even for a single row or column.
            x <- matrix(x[, , slice, drop = FALSE], nrow = dimensions[1], ncol = dimensions[2])
      } else if(!is.null(slice)){
            stop("slice is only applicable to three-dimensional arrays.")
      }
      if(is.null(mzVal)){
            mzVal <- if(length(centres) >= layer) centres[layer] else "image"
      }
      if(!(is.numeric(mzVal) || is.character(mzVal)) || length(mzVal) != 1L ||
         is.na(mzVal) || (is.numeric(mzVal) && (!is.finite(mzVal) || mzVal <= 0)) ||
         (is.character(mzVal) && !nzchar(mzVal))){
            stop("mzVal must be a single positive finite m/z or a non-empty label.")
      }
      indices <- which(!is.na(x), arr.ind = TRUE)
      if(nrow(indices) == 0L){
            stop("The selected ion image contains no non-missing pixels.")
      }
      intensity <- x[indices]
      if(any(!is.finite(intensity)) || any(intensity < 0)){
            stop("Non-missing intensities must be finite and non-negative.")
      }
      if(!any(intensity > 0)){
            stop("The selected ion image contains no positive intensities.")
      }
      coordinates <- if(orientation == "xy") indices else indices[, 2:1, drop = FALSE]
      coordinates <- data.frame(x = coordinates[, 1], y = coordinates[, 2])
      if(is.null(win)){
            win <- createSpatialWindow(coordinates, clean = FALSE, plot = FALSE, pixelSize = 1)
      } else {
            if(!inherits(win, "owin")){
                  stop("win must be a spatstat owin object.")
            }
            win <- spatstat.geom::as.polygonal(win)
            if(!all(spatstat.geom::inside.owin(coordinates$x, coordinates$y, win))){
                  stop("win must contain all non-missing pixel centres.")
            }
      }
      keep <- if(zero.rm) intensity > 0 else rep(TRUE, length(intensity))
      analytePointPattern(x = coordinates$x[keep], y = coordinates$y[keep],
                           intensity = intensity[keep], win = win, mzVals = mzVal)
}
