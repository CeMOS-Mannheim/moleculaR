#' Simulate a Spatial Point Pattern
#'
#' -- under development --
#'
#' @details
#'
#' @return
#'
#'
#' @include manualSpatstatImport.R
#'

simulatespp <- function(lambdaHotspot, lambdaBackground, 
                        intensityData,
                        radius, 
                        backgroundUpperThr = 0.75,
                        hotspotLowerThr = 0.75,
                        shape = "centralCircle",
                        backgroundLength = 100, 
                        seed = NULL,
                        gtPlot = FALSE,
                        sppPlot = TRUE){
  
  set.seed(seed)
  
  ## create owin objects - background
  bgWin <- owin(c(0, backgroundLength), c(0, backgroundLength))
  
  ## create owin objects - hotspot
  if(shape == "centralCircle"){
    if(radius > backgroundLength/2)
      stop("radius is bigger the the x-dimention of the background window. \n")
    
    if(radius > backgroundLength/2)
      stop("radius is bigger the the y-dimention of the background window. \n")
  }
  
  if(shape == "multiCircle"){
    
    xmax <- backgroundLength
    xmin <- 0
    ymax <- backgroundLength
    ymin <- 0
    
    if(any(backgroundLength/4 - radius < xmin))
      stop("one of the radii is bigger the the x-dimention of the background window. \n")
    
    if(any(backgroundLength/4*3 + radius > xmax))
      stop("one of the radii is bigger the the x-dimention of the background window. \n")
    
    if(any(backgroundLength/4 - radius < ymin))
      stop("one of the radii is bigger the the y-dimention of the background window. \n")
    
    if(any(backgroundLength/4*3 + radius > ymax))
      stop("one of the radii is bigger the the y-dimention of the background window. \n")
    
  }
  
  
  if(shape == "donut"){
    if(any(radius > backgroundLength/2))
      stop("one of the radii is bigger the the x/y-dimension of the background window. \n")

  }
  
  hotspotwin <- createHotspotWin(shape = shape, radius = radius, 
                                  backgroundLength = backgroundLength)
  

  # populate both windows with points different point intensities
  bgspp <- rpoispp(lambdaBackground, win = bgWin)
  hotspotspp <- rpoispp(lambdaHotspot, win = hotspotwin)
  
  # attach marks - background
  bgIntensities <- intensityData[intensityData <= quantile(intensityData, backgroundUpperThr)]
  bgspp$marks <- sample(bgIntensities, size = bgspp$n, replace = FALSE)
  
  # attach marks - hotspot
  hsIntensities <- intensityData[intensityData > quantile(intensityData, hotspotLowerThr)]
  hotspotspp$marks <- sample(hsIntensities, size = hotspotspp$n, replace = FALSE)
  
  
  # convert to analytePointPattern
  bgspp <- moleculaR::analytePointPattern(x = bgspp$x, y = bgspp$y, win = bgspp$window, 
                                        intensity = bgspp$marks, mzVals = "background")
  
  hotspotspp <- moleculaR::analytePointPattern(x = hotspotspp$x, y = hotspotspp$y, win = hotspotspp$window, 
                                          intensity = hotspotspp$marks, mzVals = "foreground")
  
  if(gtPlot){
    plot.owin(bgWin, 
              ylim = rev(bgWin$yrange),
              main = paste0("Simulated Groundtruth"),)
    plot.owin(hotspotwin,
              col = rgb(0,1,0,1),
              add = TRUE)  
  }
  
  
  # superimpose
  gtspp <- moleculaR::superimposeAnalytes(list(bgspp, hotspotspp))
  
  if(sppPlot){
    moleculaR::plotAnalyte(gtspp, size = 1, show.window = FALSE,
                           transpFactor = 0.5, leg.side = "right", 
                           leg.args = list(cex = 3, cex.axis = 1.25))  
  }
  
  return(gtspp)
  
}


#' Create a Hotspot Window for Simulation
#'
#'
#'
#' @details
#'
#' @return
#'
#' @export
#' @include manualSpatstatImport.R
#'
createHotspotWin <- function(shape = "centralCircle", radius = 25, backgroundLength = 100){

         hotspotwin <- switch(shape,
                              "centralCircle" = {
                                       if(length(radius) > 1)
                                                stop("for shape 'centralCircle' the radius must have 1 value. \n")

                                       disc(radius = radius, centre = c(backgroundLength/2, backgroundLength/2))
                              },
                              "multiCircle" = {

                                       if(length(radius) != 5 & length(radius) != 1)
                                                stop("for shape 'multiCircle' the radius must have either 1 or 5 values. \n")

                                       if(length(radius) == 1)
                                                radius <- rep(radius, 5)

                                       centers <- data.frame(x = c(backgroundLength/4,
                                                                   backgroundLength/4*2,
                                                                   backgroundLength/4*3,
                                                                   backgroundLength/4*3,
                                                                   backgroundLength/4),
                                                             y = c(backgroundLength/4,
                                                                   backgroundLength/4*2,
                                                                   backgroundLength/4*3,
                                                                   backgroundLength/4,
                                                                   backgroundLength/4*3))

                                       winlist <- lapply(seq(1, length(radius)), FUN = function(i){
                                                disc(radius = radius[i], centre = centers[i, ])
                                       })

                                       do.call("union.owin", winlist)

                              },

                              "donut" = {

                                       if(length(radius) != 2)
                                                stop("for shape 'donut' the radius must have either 2 values. \n")

                                       if(diff(radius) == 0)
                                                stop("for shape 'donut' the radius must have either 2 not equal values. \n")

                                       radius <- sort(radius, decreasing = TRUE)

                                       winlist <- lapply(radius, FUN = function(i){
                                                disc(radius = i, centre = c(backgroundLength/2, backgroundLength/2))
                                       })

                                       names(winlist) <- c("A", "B")

                                       do.call("setminus.owin", winlist)

                              },
                              stop("shape must be one of c('centralCircle', 'multiCircle', 'donut'). \n"))

         return(hotspotwin)


}

