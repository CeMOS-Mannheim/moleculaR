#' Analyte detection in an MSI dataset
#'
#' Detect an analyte directly in a native Cardinal MSI object. The legacy
#' 'sparseIntensityMatrix' input remains supported during migration.
#'
#' @param m the mass to be queried (in Da).
#' @param fwhm the fwhm at 'm'.
#' @param spData A Cardinal 'MSImagingArrays' or 'MSImagingExperiment' object,
#'   or a legacy 'sparseIntensityMatrix'. Cardinal input must contain known
#'   centroided spectra from a single run and no queued preprocessing.
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
#'   No detections return an empty `ppp`, matching the legacy convention.
#'
#' @details For Cardinal input, select measured peaks inclusively within three
#'   Gaussian standard deviations of `m`, where sigma = fwhm / (2 * sqrt(2 * log(2))).
#'   Gaussian weights use [gaussWeight()]. This strict interval can exclude peaks
#'   included by the legacy nearest-index boundary selection.
#'
#'   For unaligned arrays, `mean` averages measured peaks inside each spectrum's
#'   interval, including stored zeros. For shared-axis experiments, it averages
#'   all selected features including zero-filled missing peaks. These denominators
#'   can differ; arrays do not reconstruct the legacy global feature axis.
#'   Empty intervals return zero for every method. Invalid retrieved spectral
#'   values raise errors. Duplicate and unsorted measured peaks are supported.
#'
#'   All acquired pixels, including zero analyte intensities, define the spatial
#'   window. Only positive intensities become points. Multiple runs or spatial
#'   z planes must be subset explicitly with Cardinal. Call `Cardinal::process()`
#'   before extraction if preprocessing is queued. User metadata names must be
#'   unique and exclude the reserved fields `idx`, `mzVals`, and `mzConfirmed`.
#'
#' @export
#' @md
#' @include manualSpatstatImport.R
#'
searchAnalyte     = function(m, fwhm, spData, wMethod = "Gaussian", spwin = NA,
                             verifiedMasses = NA, confirmedOnly = FALSE,
                             metaData = list()) {

        if(.msiIsObject(spData)) {
              return(.searchCardinalAnalyte(m, fwhm, spData, wMethod, spwin,
                                            verifiedMasses, confirmedOnly, metaData))
        }

        if(!inherits(spData, "sparseIntensityMatrix"))
              stop("spData must be a Cardinal MSI object or legacy sparseIntensityMatrix.")




        if(length(m) > 1 | length(fwhm) > 1){
                stop("'m' length must be equal to one. \n")
        }

         if(!is.numeric(m) | !is.numeric(fwhm)){
                  stop("'m' and 'fwhm' must have numeric value. \n")
         }

         if(is.na(m) | is.na(fwhm)){
                  stop("'m' and 'fwhm' must have numeric value. \n")
         }

        if(!identical(spwin, NA)){
              owinInfo <- summary(spwin)

              if(owinInfo$npoly > 1){

                     if(any(owinInfo$nvertices <= 4)){
                            warning("The provided window 'spwin' contains several polygonal areas, ",
                            "one of which does not seem to be a tissue section with number of ",
                            "vertices less or equal to 4. \n ")
                     }

              } else{
                    if(owinInfo$nvertices <= 4){
                          warning("The provided window 'spwin' does not seem to be a tissue section ",
                                  "with number of vertices less or equal to 4. \n ")
                    }
              }

        }



       # find sigma of the gaussian --> corresponds to search window
       s             = fwhm / (2*sqrt(2*log(2))) # sigma i.e. std
       threeS         = s * 3


       #// check if listed in verifiedMasses
       if(identical(verifiedMasses, NA)) { # if identification were not supplied flag it to FALSE

               mzConfirmed      = FALSE


       } else { # identifications supplied


               # match closest
               idxConf        = MALDIquant::match.closest(m, verifiedMasses, threeS, NA_integer_)

               # record
               mzConfirmed    = ifelse(identical(idxConf, NA_integer_), FALSE, TRUE)

       }

       # window object
       if(identical(spwin, NA)){
             spwin <- as.polygonal(owin(mask = spData$coordinates))
       }

       if(confirmedOnly){ # if only confirmed detections are needed

               if(!mzConfirmed){ # and there are no confirmed detections then return empty ppp
                       return(.craeteEmptySpp(spwin))
               }
       }

       idx           <- MALDIquant::match.closest(m , spData$mzAxis, threeS, NA)


       if(!is.na(idx))
       {

              idxlwr        <- MALDIquant::match.closest((m - threeS), spData$mzAxis, threeS, idx)
              idxupr        <- MALDIquant::match.closest((m + threeS), spData$mzAxis, threeS, idx)


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

      return(ppp(x = integer(0),
                                y = integer(0),
                                marks = data.frame(idx = integer(0), intensity = numeric(0)),
                                window = spwin,
                                checkdup = FALSE,
                                drop = FALSE))
}

