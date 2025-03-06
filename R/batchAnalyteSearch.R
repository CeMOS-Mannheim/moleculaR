#' Batch Analyte detection in an MSI dataset
#'
#' This applies `moleculaR::searchAnalyte` on a given MSI dataset against a vector of provided m/z values representing
#' analytes of interest.
#'
#' @param spData an S3 object of type `sparseIntensityMatrix` holding the sparse MSI data.
#' @param fwhmObj an S3 object of type `fwhm` with the estimated fwhm data.
#' @param spwin optional, an object of type `owin`. If not given the function tries to generate the
#' spatial window out of the coordinates of all points of the dataset stored in `spData` (default behavior).
#' @param m a numeric vector of m/z values of interest.
#' @param numCores an integer, the number of cores to be used for the search, passed to `parallel::mclapply`. Not
#' supported on Windows machines.
#' @param wMethod wighting method; c("Gaussian", "sum", "max", "mean").
#' @param verifiedMasses an optional numeric vector of m/z values that are (externally) verified to
#' be real molecular entities with a certain confidence, ex. 'mz' column of a METASPACE
#' annotation result.
#' @param confirmedOnly if `TRUE`, returns detections only if confirmed by `verifiedMasses`.
#' @param verbose whether to show progress. Ignored when `numCores > 1`.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing all analyte hits identified in the
#' MSI dataset `spData` for the specified m/z values in `m` taking into account the `fwhm` information.
#'
#' @export
#'
#' @include manualSpatstatImport.R


batchAnalyteSearch <- function(spData, fwhmObj, spwin = NA, m,
                             numCores = 1L, wMethod = "Gaussian", verifiedMasses = NA,
                             confirmedOnly = FALSE, verbose = TRUE) {


      # check OS type
      if(.Platform$OS.type == "windows" & numCores > 1){
            warning("Only single-core operation is supported on windows. \n")
            numCores <- 1L
      }

      #// sort verifiedMasses
      if(!identical(verifiedMasses, NA)){

            if(!is.numeric(verifiedMasses)){
                  stop("verifiedMasses is provided but does not appear to be numeric. \n")
            }

            verifiedMasses <- sort(verifiedMasses)

      }

      #// create sp window
      if(identical(spwin, NA)){
            spwin <- createSpatialWindow(pixelCoords = spData$coordinates, plot = FALSE)
      }



      if(verbose & numCores == 1){
            pb <- utils::txtProgressBar(min = 0, max = length(m), style = 3, width = 20)
      }




      hitsList <- parallel::mclapply(X = seq(1, length(m)), mc.cores = numCores,
                                     FUN = function(i) {


                                           if(verbose & numCores == 1){
                                                 utils::setTxtProgressBar(pb, i)
                                           }


                                           # initialize empty spp objects for each adduct
                                           sppCotainer <- .craeteEmptySpp(spwin) # this function is in `searchAnalyte.R`

                                           mQuery    <- m[i]

                                           massNotNA <- !(is.na(mQuery))
                                           massInRange <- ifelse(massNotNA,
                                                                 check.in.range(mQuery, range(spData$mzAxis), FALSE),
                                                                 FALSE)

                                           if(massNotNA & massInRange) {


                                                 sppCotainer <- searchAnalyte(m = mQuery,
                                                                              fwhm = getFwhm(fwhmObj, mQuery),
                                                                              spData = spData,
                                                                              spwin = spwin,
                                                                              wMethod = wMethod,
                                                                              verifiedMasses = verifiedMasses,
                                                                              confirmedOnly = confirmedOnly,
                                                                              metaData = list())


                                           }

                                           rm(massNotNA, massInRange, mQuery)





                                           return(sppCotainer)


                                     })




      #// merge
      hitsList <- superimposeAnalytes(hitsList, spWin = spwin, check = FALSE)



      if(verbose & numCores == 1)
            close(pb)


      return(hitsList)


}




