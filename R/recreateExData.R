#' Re-create example data
#'
#' This function is for internal use only. It re-creates the example data used within the package.
#'
#' @param pathToImzml:        path to imzML file.
#' @param pathToSingleSpctr:  path to the single spectrum tsv file.
#' @param pathToSldb:         path to the swisslipids database tsv file.
#' @param pathToMtspc:        path to the Metaspace annotation csv file.
#' @param shrink:             whether to shring the msData to only retain the identified m/z in
#' the Metaspace annotation file.
#' @return
#' Saves a RData file of the example data.
#'
#' @export
#' @keywords internal
#'
recreateExData     = function(pathToImzml, pathToSingleSpctr = NULL,
                              pathToSldb = NULL, pathToMtspc = NULL,
                              shrink = TRUE) {



      #// read the file into R
      msData      = readCentrData(path = pathToImzml)


      #// single spectrum -  load a local tsv file
      if(is.null(pathToSingleSpctr)){
            pathToSingleSpctr = system.file("extdata",
                                            "pos-XIII-82492-singleSpectrum.tsv",
                                            package = "moleculaR",
                                            mustWork = TRUE)

      }
      msSpectr    = readSingleSpect(pathToSingleSpctr)


      #// load the processed swisslipids db internally
      if(is.null(pathToSldb)){
            pathToSldb = system.file("extdata", "swisslipids-speciesOnly-sep2020.tsv",
                                     package = "moleculaR", mustWork = TRUE)

      }
      sldb        = loadSwissDB(pathToSldb)


      #// load the metaspace annotations file
      if(is.null(pathToMtspc)) {
            pathToMtspc = system.file("extdata", "metaspace_annotations.csv",
                                      package = "moleculaR", mustWork = TRUE)
      }
      mtspc       = read.csv(file = pathToMtspc, skip = 2,header = TRUE, colClasses = "character")


      #// pre-processing

      #// estimate fwhm from msSpectr ----
      fwhm        = estimateFwhm(s = msSpectr, plot = FALSE)


      # bin peaks
      msData      =  MALDIquant::binPeaks(msData,
                                          tolerance = fwhmFun(400)/400, #focusing on lipids
                                          method = "relaxed")


      # filter out peaks which occur in less than 1% of the time - the
      # built-in function of MALDIquant crashes for bigger datasets
      msData      = filterPeaks(x = msData, minFreq = 0.01)


      # normalization
      #msData     = foldChangeNorm(msData) <-- this is optional


      #// create spatial window
      spwin       = spatstat::as.polygonal(spatstat::owin(mask = as.data.frame(MALDIquant::coordinates(msData))))
      spwin$bdry  = spwin$bdry[-2] # fix the small glitch

      #// shrink msData
      if(shrink){
            msData= shrinkData(x = msData, mzKeep = as.numeric(mtspc$mz), fwhmFun = fwhm$fwhmFun)

      }


      #// remove everything else
      rm(pathToImzml, pathToSingleSpctr, pathToSldb, pathToMtspc,
         shrink, msSpectr)

      #// save image to working dir
      save(list = c("msData", "sldb", "mtspc", "fwhm", "spwin"),
           file = "processed-example-Data.RData", compress = TRUE)


      cat("Done. \n")
      return(NULL)

}

