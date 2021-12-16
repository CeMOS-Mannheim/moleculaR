#' Analyte detection in an MSI dataset
#'
#' This function check if an analyte (given by m/z) is detectable in the MS data represented
#' by a 'sparseIntensityMatrix' object.
#'
#' @param m the mass to be queried (in Da).
#' @param fwhm the fwhm at 'm'.
#' @param spData an S3 object of type 'sparseIntensityMatrix' holding the sparse MSI data.
#' @param wMethod wighting method; c("Gaussian", "sum", "max", "mean").
#' @param spwin optional, an object of type `owin`. If not given the function tries to generate the
#' spatial window out of the coordinates of all points of the dataset stored in `spData` (default behavior).
#' @param verifiedMasses an optional numeric vector of m/z values that are (externally) verified to
#' be real molecular entities with a certain confidence, ex. 'mz' column of a METASPACE
#' annotation result.
#' @param confirmedOnly if `TRUE`, returns detections only if confirmed by `verifiedMasses`.
#' @param metaData optional named list with additional identifiers for the analyte
#' under study, ex. list(lipidID = "..", sumformular = "..", ..). This will be passed to the `metaData`
#' slot of the resulting `analytePointPattern` object.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing the hits for
#' `m` in `msData`. If `confirmedOnly = TRUE` then only these hits are returned.
#'
#' @export
#'
searchAnalyte     = function(m, fwhm, spData, wMethod = "Gaussian", spwin = NA,
                             verifiedMasses = NA, confirmedOnly = FALSE,
                             metaData = list()) {


        if(length(m) > 1){
                stop("'m' length must be equal to one. \n")
        }



       # find sigma of the gaussian --> corresponds to search window
       s             = fwhm / 2.355 # sigma i.e. std
       fiveS         = s * 5


       #// check if listed in verifiedMasses
       if(identical(verifiedMasses, NA)) { # if identification were not supplied flag it to FALSE

               mzConfirmed      = FALSE


       } else { # identifications supplied


               # match closest
               idxConf        = MALDIquant::match.closest(m, verifiedMasses, (s * 3), NA_integer_)

               # record
               mzConfirmed    = ifelse(identical(idxConf, NA_integer_), FALSE, TRUE)

       }

       # window object
       if(identical(spwin, NA)){
             spwin <- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = spData$coordinates))
       }

       if(confirmedOnly){ # if only confirmed detections are needed

               if(!mzConfirmed){ # and there are no confirmed detections then return empty ppp
                       return(.craeteEmptySpp(spwin))
               }
       }

       idx           <- MALDIquant::match.closest(m , spData$mzAxis, fiveS, NA)


       if(!is.na(idx))
       {

              idxlwr        <- MALDIquant::match.closest((m - fiveS), spData$mzAxis, fiveS, idx)
              idxupr        <- MALDIquant::match.closest((m + fiveS), spData$mzAxis, fiveS, idx)


              coi           <- as(spData$spmat[ , (idxlwr:idxupr), drop = FALSE], "matrix")  #columns of interest

              combinedCols  <- switch(wMethod,
                     "Gaussian" = {
                            gw <- gaussWeight(x = spData$mzAxis[(idxlwr:idxupr)],
                                          m = m,
                                          fwhm = fwhm)

                            coi           <- sweep(coi, MARGIN = 2, gw, "*")

                            rowSums(coi)

                     },
                     "sum" = {
                            rowSums(coi)
                     },
                     "mean" = {
                            rowMeans(coi)
                     },
                     "max" = {
                            apply(coi, 1, max)
                     })




              detectedIn    <- which(combinedCols > 0) # in which spectra it had a non-zero value

              if(length(detectedIn) > 0) {

                     # coordinates of these spectra
                     detectedCoord <- spData$coordinates[detectedIn, , drop = F]



                     # analytePointPattern
                     rObj <- analytePointPattern(x = detectedCoord[ , "x"],
                                                 y = detectedCoord[ , "y"],
                                                 win = spwin,
                                                 intensity = combinedCols[detectedIn],
                                                 mzVals = m,
                                                 metaData = c(list(mzConfirmed = mzConfirmed), metaData))

                     return(rObj)

              }

       }



        return(.craeteEmptySpp(spwin)) # return empty ppp object







}

.craeteEmptySpp <- function(spwin){

      # creates an empty spp compatible with the internal representation
      # of moleculaR

      return(spatstat.geom::ppp(x = integer(0),
                                y = integer(0),
                                marks = data.frame(idx = integer(0), intensity = numeric(0)),
                                window = spwin,
                                checkdup = FALSE,
                                drop = FALSE))
}

