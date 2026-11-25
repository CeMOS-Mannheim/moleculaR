#' Batch Analyte detection in an MSI dataset
#'
#' This applies `moleculaR::searchAnalyte` on a given MSI dataset against a vector of provided m/z values representing
#' analytes of interest.
#'
#' @param spData A native Cardinal `MSImagingArrays` or `MSImagingExperiment`,
#'   or a temporary legacy `sparseIntensityMatrix`.
#' @param fwhmObj an S3 object of type `fwhm` with the estimated fwhm data.
#' @param spwin optional, an object of type `owin`. If not given the function tries to generate the
#' spatial window out of the coordinates of all points of the dataset stored in `spData` (default behavior).
#' @param m a numeric vector of m/z values of interest.
#' @param numCores Number of workers. Cardinal uses serial execution by default
#'   and a cross-platform BiocParallel SnowParam for values greater than one.
#'   The legacy path uses `parallel::mclapply` and is serial on Windows.
#' @param wMethod wighting method; c("Gaussian", "sum", "max", "mean").
#' @param verifiedMasses an optional numeric vector of m/z values that are (externally) verified to
#' be real molecular entities with a certain confidence, ex. 'mz' column of a METASPACE
#' annotation result.
#' @param confirmedOnly if `TRUE`, returns detections only if confirmed by `verifiedMasses`.
#' @param verbose Logical; show pixel-block progress for Cardinal input. The
#'   legacy path ignores this when `numCores > 1`.
#' @param chunkSize Cardinal input only: positive integer number of pixels per
#'   block. NULL uses Cardinal's configured chunk size or number of chunks.
#' @param BPPARAM Cardinal input only: optional BiocParallel backend. Leave
#'   numCores at one when supplying this argument. Caller-started backends stay
#'   running; backends started by this call are stopped on exit.
#' @param cache Cardinal input only: optional unlocked environment storing the
#'   last sparse extraction. Reuse requires identical data, processing, run/ROI,
#'   masses, widths, method and backing-file statistics. Keep backing files
#'   immutable; clear the environment after external edits preserving timestamps.
#'
#' @details Cardinal spectra are read once per batch in pixel blocks. Extraction
#'   follows [searchAnalyte()] including strict inclusive mass windows and the
#'   representation-specific mean denominator. Input order and duplicate masses
#'   are preserved in detected-analyte metadata. NA targets are skipped; other
#'   invalid masses and invalid FWHM predictions raise errors. Out-of-range
#'   centres are evaluated when their windows overlap measured peaks. Empty or
#'   undetected batches return an empty ppp. Output memory grows with detections;
#'   no complete spectral or dense pixel-by-query matrix is constructed.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing all analyte hits identified in the
#' MSI dataset `spData` for the specified m/z values in `m` taking into account the `fwhm` information.
#'
#' @export
#'


batchAnalyteSearch <- function(spData, fwhmObj, spwin = NA, m,
                             numCores = 1L, wMethod = "Gaussian", verifiedMasses = NA,
                             confirmedOnly = FALSE, verbose = TRUE,
                             chunkSize = NULL, BPPARAM = NULL, cache = NULL) {

      if (.msiIsObject(spData)) {
            return(.searchCardinalAnalytes(spData, fwhmObj, spwin, m, numCores,
                   wMethod, verifiedMasses, confirmedOnly, verbose, chunkSize,
                   BPPARAM, cache))
      }
      if (!is.null(chunkSize) || !is.null(BPPARAM) || !is.null(cache)) {
            stop("chunkSize, BPPARAM and cache require Cardinal input.")
      }


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
                                                                 spatstat.utils::check.in.range(mQuery, range(spData$mzAxis), FALSE),
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




