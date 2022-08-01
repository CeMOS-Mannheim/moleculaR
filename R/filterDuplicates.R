#' Filtersout duplicated analytes
#'
#' Filters out duplicated analytes based on their m/z value. Normally used for the case
#' of CPPMs i.e. `spp` objects carrying more than one analyte.
#'
#' @param spp: 	   a spatial point pattern of type 'ppp' with its `marks` being a dataframe with
#' two columns `idx`and `intensity`.
#'
#' @return
#'
#' @export
#'

filterDuplicates <- function(spp){

         torm <- duplicated(spp$metaData$mzVals)
         tokp <- !torm


         if(length(which(torm)) == 0){
                  return(spp)
         }

         kpidx <- spp$metaData$idx[which(tokp)]

         spp <- subsetAnalytes(spp, idx %in% kpidx)

         return(spp)

}
