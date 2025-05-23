#' Re-create example data
#'
#' This function is for internal use only. It re-creates the example data used within the package.
#'
#' @param pathToImzml:        path to imzML file.
#' @param pathToSingleSpctr:  path to the single spectrum tsv file.
#' @param pathToSldb:         path to the swisslipids database tsv file.
#' @param pathToMtspc:        path to the Metaspace annotation csv file.
#' @param shrink:             whether to shrink the msData to only retain the identified m/z in
#' the Metaspace annotation file.
#' @param saveTo:             where to save the resulting example data.
#' @return
#' Saves a RData file of the example data.
#'
#' @export
#' @keywords internal
#'
recreateExData     <- function(pathToImzml, pathToSingleSpctr = NULL,
                              pathToSldb = NULL, pathToMtspc = NULL,
                              shrink = TRUE, saveTo = getwd()) {



      #// read the file into R
      cat("loading imzML data .. \n")
      msData      <- readCentrData(path = pathToImzml)



      #// single spectrum -  load a local tsv file
      if(is.null(pathToSingleSpctr)){
            pathToSingleSpctr <- system.file("extdata",
                                            "pos-XIII-82492-singleSpectrum.tsv",
                                            package = "moleculaR",
                                            mustWork = TRUE)

      }
      cat("reading single spectrum  .. \n")
      msSpectr    <- readSingleSpect(pathToSingleSpctr)


      #// load the processed swisslipids db internally
      if(is.null(pathToSldb)){
            pathToSldb <- system.file("extdata", "swisslipids-speciesOnly-sep2020.tsv",
                                     package = "moleculaR", mustWork = TRUE)

      }
      cat("loading SwissLipids data .. \n")
      sldb        <- loadSwissDB(pathToSldb)


      #// load the metaspace annotations file
      if(is.null(pathToMtspc)) {
            pathToMtspc <- system.file("extdata", "metaspace_annotations.csv",
                                      package = "moleculaR", mustWork = TRUE)
      }
      mtspc       <- read.csv(file = pathToMtspc, skip = 2,header = TRUE, colClasses = "character")


      #// pre-processing

      #// estimate fwhm from msSpectr ----
      cat("preprocessing: fwhm estimation .. \n")
      fwhmObj        <- estimateFwhm(s = msSpectr, plot = FALSE)


      # bin peaks
      cat("preprocessing: peak binning .. \n")
      msData      <-  MALDIquant::binPeaks(msData,
                                          tolerance = getFwhm(fwhmObj, 400)/400, #focusing on lipids
                                          method = "relaxed")


      # filter out peaks which occur in less than 1% of the time - the
      # built-in function of MALDIquant crashes for bigger datasets
      cat("preprocessing: peak filtration .. \n")
      msData      <- filterPeaks(x = msData, minFreq = 0.01)


      # normalization
      #msData     <- foldChangeNorm(msData) <-- this is optional


      #// create spatial window
      #spwin       <- as.polygonal(owin(mask = as.data.frame(MALDIquant::coordinates(msData))))
      #spwin$bdry  <- spwin$bdry[-2] # fix the small glitch

      #// shrink msData
      if(shrink){
            cat("preprocessing: shrinking data .. \n")
            msData<- .shrinkData(x = msData, mzKeep = as.numeric(mtspc$mz), fwhmObj = fwhmObj)

      }

      #// create sparseIntensityMatrix representation
      cat("preprocessing: creating sparseIntensityMatrix .. \n")
      spData      <- createSparseMat(msData)


      #// remove everything else
      rm(pathToImzml, pathToSingleSpctr, pathToSldb, pathToMtspc,
         shrink, msSpectr)

      #// save image to working dir
      save(list = c("msData", "sldb", "mtspc", "fwhmObj"),
           file = file.path(saveTo, "processed-example-Data.RData"), compress = TRUE)


      cat("Done. \n")
      return(NULL)

}


#' Shrink data
#'
#' This function is for internal use only. It tries to shrink MSI data (in terms of size) for including it within
#' a package.
#'
#' @param x: 	   Dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param mzKeep: numeric vector, m/z values to retain.
#' @param fwhmObj: the 'fwhm' S3 object. See \code{?moleculaR::estimateFwhm}.
#' @return
#' Filtered list of \code{MALDIquant::MassPeaks} objects.
#'
#' @export
#' @keywords internal
#'
.shrinkData     <- function(x, mzKeep, fwhmObj) {

   if(!MALDIquant::isMassPeaksList(x)) {
      stop("x must be a list of MassPeaks objects.\n")
   }

   mzKeep      <- sort(mzKeep)
   x           <- lapply(x, FUN = function(i) {

      xm       <- MALDIquant::mass(i)
      xkpIdx   <- MALDIquant::match.closest(x = mzKeep,
                                           table = xm,
                                           tolerance = getFwhm(fwhmObj, mzKeep))
      xkpIdx   <- xkpIdx[!is.na(xkpIdx)]

      mt       <- list()
      mt$imaging$pos <- MALDIquant::metaData(i)$imaging$pos # keep only pixel coordinates


      MALDIquant::createMassPeaks(mass = xm[xkpIdx],
                                  intensity = MALDIquant::intensity(i)[xkpIdx],
                                  snr = MALDIquant::snr(i)[xkpIdx],
                                  metaData = mt)

   })


   return(x)

}

