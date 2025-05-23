% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSpatialWindow.R
\name{createSpatialWindow}
\alias{createSpatialWindow}
\title{Create Spatial Window of MSI Data}
\usage{
createSpatialWindow(pixelCoords, isPolygon = FALSE, clean = FALSE, plot = TRUE)
}
\arguments{
\item{pixelCoords:}{a data frame with x and y columns containing the coordinates
of all pixels within a given MSI data. A `matrix` will be coerced into a
data frame.}

\item{isPolygon:}{a logical, whether `pixelCoords` represents the coordinates
of a polygon or not.}

\item{clean:}{logical, whether to fill holes created by missing pixels.}

\item{plot:}{whether to plot.}
}
\value{
A spatial window object of type `spatstat.geom::owin`.
}
\description{
This function creates `spatstat.geom::owin` polygonal spatial window which represents
the boundaries of the tissue under study.
}
