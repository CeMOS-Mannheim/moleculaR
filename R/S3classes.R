## Classes
##
## Under the hood, moleculaR relies on MALDIquant, MALDIquantForeign,
## Matrix and spatstat packages.It inherits classes from 'Matrix' and
## 'spatstat' packages. For more info please check the documentation
## the respective packages.
##
##

# to document this -> https://r-pkgs.org/man.html#man-functions


#' fwhm Class Constructor
#'
#' This constructs a 'fwhm' object which stores information about the calculated full width half maximum.
#'
#' @param fwhmInterpolator A function that is used to find fwhm value for any given m/z.
#' @param  peaks  The m/z values of the peaks of a single spectrum.
#' @param  fwhmVals   The corresponding fwhm values in Da.
#'
#' @return
#' An S3 object 'fwhm'.
#'
#' @export
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

print.fwhm <- function(obj) {

      cat("S3 object 'fwhm' holding fwhm data of a certain mass spectrum. \n ")
      cat("fwhm data is based on ", length(obj$peaks), "fwhm values,.\n")

}

getFwhm <- function(obj, mzVals) {
      UseMethod("getFwhm")
}

getFwhm.fwhm <- function(obj, mzVals){
      obj$fwhmInterpolator(mzVals)
}



#' sparseIntensityMatrix Class Constructor
#'
#' An S3 class to hold mass spectrometry imaging (MSI) intensity data in a sparse representation ('Matrix::dgCMatrix' class)
#' with few additional slots. It is used to speed-up computation.
#'
#' @param spmat A sparse matrix of type 'Matrix::dgCMatrix' holding MSI intensity data.
#' @param mzAxis A numeric, the overall m/z axis for the experiment.
#' @param coordinates  A two-column matrix, the coordinates of the spectra which are stored
#' in rows of 'spmat'.
#'
#' @return
#' An S3 object 'sparseIntensityMatrix'.
#'
#' @export
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

if(class(obj$coordinates) != "matrix"){
      stop("slot 'coordinates' must be a matrix. \n")
}

if(ncol(obj$coordinates) != 2){
      stop("slot 'coordinates' must be a matrix of two columns. \n")
}

if(nrow(obj$coordinates) != obj$spmat@Dim[1]){
      stop("number of rows of @coordinates must be equal to the number of rows of the sparse matrix.\n")
}


class(obj) <- "sparseIntensityMatrix"

return(obj)

}

print.sparseIntensityMatrix <- function(obj) {

      cat("S3 object 'sparseIntensityMatrix' holding sparse intensity data of an MSI dataset. \n ")
      cat("spmat: sparse matrix of type 'Matrix::dgCMatrix' with ", nrow(obj$spmat),
          "rows (spectra/pixels) and ", ncol(obj$spmat), "columns (m/z bins) .\n")
      cat("mzAxis: a numeric vector holding the m/z axis = number of columns of spmat.\n")
      cat("coordinates: a 2-column matrix holding the spectra/pixel coordinates. \n")

}




#' analytePointPattern Class Constructor
#'
#' This is a class to extend 'spatstat::ppp' class only enforcing specific point-marks. It also
#' adds an 'mzVals' slot to declare which analytes are constructing a given 'analytePointPattern' object.
#' It represents the spatial point pattern of an analyte (or set of
#' analytes) with marks being a data.frame of one column c("intensity").
#'
#' @param spp a spatial point pattern of type 'spatstat::ppp'.
#' @param x,y,intensity,win optional, either these arguments or 'spp' have to be supplied.
#' @param mzVals a numeric vector holding the m/z values that represent the analyte under study.
#'
#' @return
#' An S3 object of type 'ppp' and 'analytePointPattern'.
#'
#' @export
#'
#'

analytePointPattern <- function(spp = NA, x = NA, y = NA, intensity = NA, win = NA, mzVals) {


      # integrity checks
      if(all(is.na(spp), any(is.na(x, y, intensity, win)))){
         stop("either 'spp' or 'x,y,intensity,win' must be supplied. \n")
      }

      if(length(mzVals) < 1){
            stop("'mzVals' slot length must be more than one.\n")
      }

      if(class(spp) != "ppp"){
         stop("slot 'spp' must be of type 'spatstat::ppp' \n")
      }

      if(class(spp$marks) != "data.frame"){
         stop("marks of 'spp' must contain a data.frame of single column named 'intensity'. \n")
      }

      if(identical(colnames(spp$marks), "intensity")){
         stop("marks of 'spp' must contain a data.frame of single column named 'intensity'. \n")
      }


      # definition
      if(is.na(spp)){


         spp <- spatstat::ppp(x = x, y = y,
                              marks = data.frame(intensity = intensity),
                              window = win, checkdup = FALSE, drop = FALSE)


      }

      spp$mzVals <- mzVals
      class(spp) <- c(class(spp), "analytePointPattern")

      return(spp)

}






#' molProbMap Class Constructor
#'
#' This is a class to represent the calculated molecular probability maps.
#'
#' @param bw   A numeric, the calculated Gaussian bandwidth
#' @param sppMoi A 'analytePointPattern' object representing a given analyte (metabolite-of-interest).
#' @param csrMoi A 'analytePointPattern' object representing the corresponding CSR model.
#' @param rhoMoi An 'spatstat::im' object storing the density image of 'sppMoi' after applying KDE.
#' @param rhoCsr An 'spatstat::im' object storing the density image of 'csrMoi' after applying KDE.
#' @param hotspotpp A 'analytePointPattern' object holding only the points that are located within the hotspot.
#' @param hotspotIm An 'spatstat::im' object storing the corresponding pixellated image of the 'hotspotpp' slot.
#' @param hotspotMask An 'spatstat::owin' object storing the corresponding window of the 'hotspotpp' slot.
#' @param coldspotpp A 'analytePointPattern' object holding only the points that are located within the coldspot.
#' @param coldspotIm An 'spatstat::im' object storing the corresponding pixellated image of the 'coldspotpp' slot.
#' @param coldspotMask An 'spatstat::owin' object storing the corresponding window of the 'coldspotpp' slot.
#'
#' @return
#' An S3 object of type 'molProbMap'.
#'
#' @export
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

print.molProbMap <- function(obj) {

      cat("S3 object 'molProbMap' holding molecular probability maps of a given analyte in an MSI dataset. \n ")
      cat("spp: spatial point pattern of type 'spatstat::ppp' with ", obj$spp$n, " points/events.\n")
      cat("mzVals: m/z values that represent the analyte under study = ", round(object$mzVals, 4), " .\n")

}


#' lipidSearchList Class Constructor
#'
#' A class to re-arrange lipid listings in the SwissLipid database to facilitate computational
#' serach. It splits the full database into lipid class-oriented lists.
#'
#' @param lipidList a data.frame which lists entries of lipids.
#' @param hitsList a data.frame holding the lipids that were detected.
#' @param allClasses a character vector summarizing all classes taking part in the lipid search.
#'
#' @return
#' An S3 object of type 'lipidSearchList' .
#'
#' @export
#'
#'

lipidSearchList <- function(lipidList, hitsList, allClasses) {


   # definition
   l <- list(lipidList = lipidList, hitsList = hitsList, allClasses = allClasses)
   class(l) <- "lipidSearchList"

   return(l)

}
