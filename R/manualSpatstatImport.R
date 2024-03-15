#' import spatstat sub-packages
#'
#' This is done because spatstat is in constant development and on
#' several occasions some functions jumped from one sub-package to
#' another (there were also newly created sub-packages). To avoid attaching the entire
#' namespace of spatstat and its sub-packages (by listing it under "Depends"),
#' the below loader script is called when the package is loaded
#' to make it possible to work with tha package without attaching it
#' via "library()". Note that this adds the requirement of listing
#' all spatstat functions needed for normal functioning of moleculaR in the
#' function below, but will (hopefully) ensure moleculaR's compatibility with
#' all spatstat versions.
#'
#' This is a temporary solution until spatstat reaches a steady state.
#'
#' @name manualSpatstatImport
#'



# check which spatstat packages are installed
installedPckgs <- rownames(installed.packages())
spatstatPckgs <- grep("spatstat", installedPckgs, value = TRUE)

# enumerate the functions used from spatstat
spatstatFuncs <- c("ppp", "as.ppp", "is.ppp", "plot.ppp",
                   "owin", "as.owin", "area.owin", "plot.owin",
                   "as.polygonal", "as.mask",
                   "rotate", "rpoint", "npoints",
                   "density.ppp", "spatdim", "check.1.integer",
                   "coords", "colourmap", "to.transparent", "symbolmap",
                   "im", "plot.im", "pixellate", "eval.im", "as.data.frame.im", "as.im.data.frame",
                   "subset.ppp", "superimpose", "superimpose.ppp",
                   "is.empty.ppp", "intersect.owin", "rpoispp", "disc", "as.im",
                   "erosion", "Window", "Smooth.im", "Area.xypolygon", "check.in.range",
                   "blur", "resolve.1.default", "safeDevCapabilities",
                   "owinInternalRect")

correspondPckgs <- setNames(vector("character", length(spatstatFuncs)), spatstatFuncs)

# functions belong to which spatstat sub-package?
for(ifunc in spatstatFuncs){
      for(ipckg in spatstatPckgs){
            allFuncs <- getNamespaceExports(ipckg)
            if(ifunc %in% allFuncs){
                  correspondPckgs[ifunc] <- ipckg
            }
      }
}

# for older spatstat versions some internal functions above might not be available
correspondPckgs <- correspondPckgs[correspondPckgs != ""]

for(ifunc in 1:length(correspondPckgs)){
      # assign(names(correspondPckgs)[ifunc],
      #        eval(parse(text = paste0(correspondPckgs[ifunc], "::", names(correspondPckgs)[ifunc]))))
      import::here(.from=correspondPckgs[ifunc],
                   names(correspondPckgs)[ifunc],
                   .character_only = TRUE)
}









