## Classes
##
## Under the hood, moleculaR relies on MALDIquant, MALDIquantForeign,
## Matrix and spatstat packages.It inherits classes from 'Matrix' and
## 'spatstat' packages. For more info please check the documentation
## the respective packages.
##
##
## to document this -> https://r-pkgs.org/man.html#man-functions

#' @include manualSpatstatImport.R
NULL

#' fwhm Class Constructor
#'
#' This constructs a 'fwhm' object which stores information about the calculated full width half maximum.
#' This function is for internal use and normally should not be called by the user.
#'
#' @param fwhmInterpolator A function that is used to find fwhm value for any given m/z.
#' @param  peaks  The m/z values of the peaks of a single spectrum.
#' @param  fwhmVals   The corresponding fwhm values in Da.
#'
#' @return
#' An S3 object 'fwhm' with the following entries \cr
#' \itemize{
#'   \item fwhmInterpolator: A function that is used to find fwhm value for any given m/z.
#'   \item peaks: The m/z values of the peaks of a single spectrum that were used for fwhm calculation.
#'   \item fwhmVals: The corresponding fwhm values in Da.
#' }
#'
#' @export
#' @keywords internal
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


#' @export
print.fwhm <- function(obj) {

      cat("S3 object 'fwhm' holding fwhm data of a certain mass spectrum. \n ")
      cat("fwhm data is based on ", length(obj$peaks), "fwhm values,.\n")

}

#' @export
plot.fwhm <- function(obj, ...) {

   r = range(obj$peaks)
   qp = seq(r[1], r[2])
   plot(x = obj$peaks, y = obj$fwhmVals,
        main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm", ...)
   lines(x = qp, y = getFwhm(obj, qp), col = "green", lwd = 2)

}

#' Get fwhm for a given m/z value
#'
#' This is a method for S3 class `fwhm` to interpolate fwhm values.
#'
#' @param obj S3 object of type `fwhm`.
#' @param mzVals m/z value (Da) at which to estimate the fwhm.
#'
#' @return returns a numeric, estimated fwhm at the given m/z.
#'
#' @export
getFwhm <- function(obj, mzVals) {
   if(class(obj) != "fwhm")
      stop("obj must be an S3 object of type 'fwhm'\n")

   obj$fwhmInterpolator(mzVals)
}





#' sparseIntensityMatrix Class Constructor
#'
#' An S3 class to hold mass spectrometry imaging (MSI) intensity data in a sparse representation ('Matrix::dgCMatrix' class)
#' with few additional slots. It is used to speed-up computation. This constructor is for
#' internal use and normally should not be called by the user.
#'
#' @param spmat A sparse matrix of type 'Matrix::dgCMatrix' holding MSI intensity data.
#' @param mzAxis A numeric, the overall m/z axis for the experiment.
#' @param coordinates  A two-column data.frame, the coordinates of the spectra which are stored
#' in rows of 'spmat'.
#'
#' @return
#' An S3 object 'sparseIntensityMatrix' with the following entries \cr
#' \itemize{
#'   \item spmat: A sparse matrix of type 'Matrix::dgCMatrix' holding MSI intensity data.
#'   \item mzAxis: A numeric, the overall m/z axis for the experiment.
#'   \item coordinates: A two-column data.frame, the coordinates of the spectra which are stored in rows of 'spmat'.
#' }
#'
#' @export
#' @keywords internal
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

if(class(obj$coordinates) != "data.frame"){
      stop("slot 'coordinates' must be a data.frame. \n")
}

if(ncol(obj$coordinates) != 2){
      stop("slot 'coordinates' must be a data.frame of two columns. \n")
}

if(nrow(obj$coordinates) != obj$spmat@Dim[1]){
      stop("number of rows of @coordinates must be equal to the number of rows of the sparse matrix.\n")
}


class(obj) <- "sparseIntensityMatrix"

return(obj)

}

print.sparseIntensityMatrix <- function(obj) {

      cat("S3 object 'sparseIntensityMatrix' holding sparse intensity data of an MSI dataset. \n")
      cat("spmat: sparse matrix of type 'Matrix::dgCMatrix' with ", nrow(obj$spmat),
          "rows (spectra/pixels) and ", ncol(obj$spmat), "columns (m/z bins) .\n")
      cat("mzAxis: a numeric vector holding the m/z axis = number of columns of spmat.\n")
      cat("coordinates: a 2-column data.frame holding the spectra/pixel coordinates. \n")

}




#' analytePointPattern Class Constructor
#'
#' This is a class to extend 'ppp' class only enforcing specific point-marks. It also
#' adds a 'metaData' slot to declare which analytes are constructing a given 'analytePointPattern' object.
#' It represents the spatial point pattern of an analyte (or set of
#' analytes) with marks being a data.frame of two columns c('idx', 'intensity'). This constructor is for
#' internal use and normally should not be called by the user.
#'
#' @param spp a spatial point pattern of type 'ppp' with its `marks` being a dataframe with
#' two columns `idx`and `intensity`.
#' @param x,y,intensity,win optional, either these arguments or 'spp' have to be supplied.
#' @param metaData a named list having additional attributes describing the individual masses in `mzVals`.
#'
#' @details `mzVals` and `metaData` are combined into a data.frame which becomes a class slot `metaData`.
#' This contains the `mzVals` and possible other meta data with each row having a unique identifier `idx`
#' which then link `metaData` with `spp$marks`.
#'
#' @return
#' An S3 object of type 'ppp' and 'analytePointPattern' with the following entries \cr
#' \itemize{
#'   \item regular `ppp` entries: see `spatstat` documentation for `ppp` class.
#'   \item metaData: A data frame with meta data describing the points in `spp$marks`.
#' }
#'
#' @export
#' @keywords internal
#'
#'

analytePointPattern <- function(spp = NA, x = NA, y = NA, win = NA, intensity = NA, mzVals, metaData = list()) {

      #>>> to do: move mzVals to metaData > affected: searchAnalyte, superimposeAnalytes, ...

      # integrity checks
      if(all(identical(spp, NA), any(is.na(x), is.na(y), is.na(intensity), is.na(win)))){
         stop("either 'spp' or 'x,y,intensity,win' must be supplied. \n")
      }

      if(!is.list(metaData)){
         stop("metaData must be a named list. \n")
      }

      if(length(mzVals) < 1){
         stop("mzVals must be a numeric vector with at least one element. \n")
      }


      # definition
      idx <- as.numeric(Sys.time()) - prod(sample(seq(1,99), 3)) # primary key; id to connect to metaData table with marks (saves memory)

      if(identical(spp, NA)){


         spp <- ppp(x = x, y = y,
                              marks = data.frame(idx = rep(idx, length(intensity)), intensity = intensity),
                              window = win, checkdup = FALSE, drop = FALSE)


      } else {

         if(class(spp) != "ppp"){
            stop("slot 'spp' must be of type 'ppp' \n")
         }

         if(class(spp$marks) != "data.frame"){
            stop("marks of 'spp' must contain a data.frame of at least one column named 'intensity'.\n")
         }

         if(!("intensity" %in% colnames(spp$marks))){
            stop("marks of 'spp' must contain a data.frame of at least one column named 'intensity'.\n")
         }

         spp$marks <- data.frame(idx = rep(idx, length(spp$marks$intensity)), intensity = spp$marks$intensity)

      }

      # remove rejected points
      spp <- as.ppp(spp)

      # metaData object
      #mtdt <- data.frame(idx = idx, mzVals = mzVals, stringsAsFactors = FALSE)
      mtdt <- list(idx = idx, mzVals = mzVals)

      #if(length(metaData) > 0 & is.list(metaData)){
         #metaData <- metaData[!sapply(metaData, is.null)] # to remove residual null entries
         spp$metaData <- as.data.frame(c(mtdt, metaData),
                                       stringsAsFactors = FALSE)
      # } else {
      #    spp$metaData <- mtdt
      # }

      class(spp) <- c(class(spp), "analytePointPattern")

      return(spp)

}









#' molProbMap Class Constructor
#'
#' This is an S3 class to represent the calculated molecular probability maps. This
#' class is only used with objects returned by `porbMap` function and should not be
#' called explicitly by the user.
#'
#' @param bw   A numeric, the calculated Gaussian bandwidth
#' @param sppMoi A 'analytePointPattern' object representing a given analyte (metabolite-of-interest).
#' @param csrMoi A 'analytePointPattern' object representing the corresponding CSR model.
#' @param rhoMoi An 'im' object storing the density image of 'sppMoi' after applying KDE.
#' @param rhoCsr An 'im' object storing the density image of 'csrMoi' after applying KDE.
#' @param hotspotpp A 'analytePointPattern' object holding only the points that are located within the hotspot.
#' @param hotspotIm An 'im' object storing the corresponding pixellated image of the 'hotspotpp' slot.
#' @param hotspotMask An 'owin' object storing the corresponding window of the 'hotspotpp' slot.
#' @param coldspotpp A 'analytePointPattern' object holding only the points that are located within the coldspot.
#' @param coldspotIm An 'im' object storing the corresponding pixellated image of the 'coldspotpp' slot.
#' @param coldspotMask An 'owin' object storing the corresponding window of the 'coldspotpp' slot.
#'
#' @return
#' An S3 object of type 'molProbMap'.
#'
#' @export
#' @keywords internal
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

#' @export
print.molProbMap <- function(obj) {

      cat("S3 object 'molProbMap' holding molecular probability maps of a given analyte in an MSI dataset. \n ")

}





#' lipidSearchList Class Constructor
#'
#' A class to re-arrange lipid listings in the SwissLipid database to facilitate computational
#' serach. It splits the full database into lipid class-oriented lists. This function is for
#' internal use and normally should not be called by the user.
#'
#' @param lipidList a data.frame which lists entries of lipids.
#' @param hitsList a data.frame holding the lipids that were detected.
#' @param allClasses a character vector summarizing all classes taking part in the lipid search.
#'
#' @return
#' An S3 object of type 'lipidSearchList'.
#'
#' @export
#' @keywords internal
#'

lipidSearchList <- function(lipidList, hitsList, allClasses) {


   # definition
   l <- list(lipidList = lipidList, hitsList = hitsList, allClasses = allClasses)
   class(l) <- "lipidSearchList"

   return(l)

}
