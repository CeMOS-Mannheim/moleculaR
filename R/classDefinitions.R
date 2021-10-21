# ##__ S4 class definitions for future versions, see S3classes.R instead __##
#
# ## issue: spatstat does not export ppp, owin and im classes (these are S3 classes).
# ## solutions: check the most recent versions.
#
#
# ## Classes
# ##
# ## Under the hood, moleculaR relies on MALDIquant, MALDIquantForeign,
# ## Matrix and spatstat packages.It inherits classes from 'Matrix' and
# ## 'spatstat' packages. For more info please check the documentation
# ## the respective packages.
# ##
# ##
#
# # to document this -> https://r-pkgs.org/man.html#man-functions
#
#
# # fwhmData Class
# #
# # This is a class that stores information about the calculated full width half maximum.
# #
# # @slot fwhmInterpolator A function that is used to find fwhm value for any given m/z.
# # @slot peaks  The m/z values of the peaks of a single spectrum.
# # @slot fwhmVals   The corresponding fwhm values in Da.
# #
# # @export
#
# setClass("fwhmData",
#          slots = list(fwhmInterpolator = "function", peaks = "numeric", fwhmVals = "numeric"),
#          prototype = list(peaks = NA_real_, fwhmVals = NA_real_))
#
# setValidity("fwhmData", function(object){
#
#
#       if(length(object@peaks) != length(object@fwhmVals)){
#             return("lengths of peaks and fwhmVals slots have to be equal. \n")
#       }
#
#
#       return(TRUE)
# })
#
# setGeneric(name="fwhm",
#            def=function(x, mzVals)
#            {
#                  standardGeneric("fwhm")
#            }
# )
#
# setMethod(f = "fwhm",
#           signature = "fwhmData",
#           definition = function(x, mzVals)
#           {
#                 x@fwhmInterpolator(mzVals)
#           }
# )
#
#
#
# # sparseIntensityMatrix Class
# #
# # This is a class to extend 'Matrix::dgCMatrix' class with few additional slots. It
# # stores MSI spectral intensities and is used to speed-up computation.
# #
# # @slot mzAxis A numeric, the overall m/z axis for the experiment.
# # @slot coordinates  A two-column matrix, the coordinates of the spectra (which are stored
# # in rows) of the MSI experiment.
# #
# # @importClassesFrom  Matrix dgCMatrix
# # @export
# #
# setClass("sparseIntensityMatrix",
#          slots = list(mzAxis = "numeric", coordinates = "matrix"),
#          contains = "dgCMatrix")
#
# setValidity("sparseIntensityMatrix", function(object){
#
#      if(length(object@mzAxis) != object@Dim[2]){
#         return("@mzAxis slot length must be equal to the number of columns of the sparse data matrix.\n")
#      }
#
#      if(class(object@coordinates) != "matrix"){
#          return("slot 'coordinates' must be a matrix. \n")
#      }
#
#       if(ncol(object@coordinates) != 2){
#          return("slot @coordinates must be a matrix of two columns. \n")
#       }
#
#      if(nrow(object@coordinates) != object@Dim[1]){
#             return("number of rows of @coordinates must be equal to the number of rows of the sparse matrix.\n")
#       }
#
#       return(TRUE)
# })
#
# require("spatstat", quietly = TRUE)
#
# # analytePointPattern Class
# #
# # This is a class to extend 'spatstat::ppp' class only enforcing specific point-marks. It also
# # adds an 'mzVals' slot to declare which analytes are constructing a given 'analytePointPattern' object.
# # It represents the spatial point pattern of an analyte (or set of
# # analytes) with marks being a data.frame of one column c("intensity").
# #
# # @import spatstat
# # @export
# #
# # @importClassesFrom spatstat ppp
# #
#
# setClass("analytePointPattern",
#          slots = list(mzVals = "numeric"),
#          contains = "ppp")
#
#
#
# # https://stackoverflow.com/questions/46851688/in-r-how-can-i-subclass-an-object-from-an-imported-class-in-my-package
# # https://stackoverflow.com/questions/18544006/how-do-i-indicate-collate-order-in-roxygen2
#
# setValidity("analytePointPattern", function(object){
#
#       if(class(object@marks) != "data.frame"){
#         return("marks must contain a data.frame of single column named 'intensity'. \n")
#       }
#
#       if(identical(colnames(object@marks), "intensity")){
#         return("marks must contain a data.frame of single column named 'intensity'. \n")
#       }
#
#       return(TRUE)
# })
#
#
# # molProbMap Class
# #
# # This is a class to represent molecular probability maps.
# #
# # @slot bw   A numeric, the calculated Gaussian bandwidth
# # @slot sppMoi A 'analytePointPattern' object representing a given analyte (metabolite-of-interest).
# # @slot csrMoi A 'analytePointPattern' object representing the corresponding CSR model.
# # @slot rhoMoi An 'spatstat::im' object storing the density image of 'sppMoi' after applying KDE.
# # @slot rhoCsr An 'spatstat::im' object storing the density image of 'csrMoi' after applying KDE.
# # @slot hotspotpp A 'analytePointPattern' object holding only the points that are located within the hotspot.
# # @slot hotspotIm An 'spatstat::im' object storing the corresponding pixellated image of the 'hotspotpp' slot.
# # @slot hotspotMask An 'spatstat::owin' object storing the corresponding window of the 'hotspotpp' slot.
# #
# #
# #
# # @import spatstat
# # @export
# #
# setClass("molProbMap",
#          slots = list(bw = "numeric",
#                       sppMoi = "analytePointPattern",
#                       csrMoi = "analytePointPattern",
#                       rhoMoi = "im",
#                       rhoCsr = "im",
#                       hotspotpp = "analytePointPattern",
#                       hotspotIm = "im",
#                       hotspotMask = "owin",
#                       coldspotpp = "analytePointPattern",
#                       coldspotIm = "im",
#                       coldspotMask = "owin"))
#
