#' Superimpose analytePointPatterns
#'
#' Superimpose any number of point patterns of type `ppp` and `analytePointPatterns`.
#'
#' @param pppObjs  a list of any number of `ppp` and `analytePointPatterns` objects.
#' @param spWin the spatial window where all points will be projected.
#' @param check logical, passed to `spatstat::superimpose`.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing all supplied analytes.
#'
#' @export
#'

superImposeAnalytes <- function(pppObjs, spWin = NULL, check = TRUE){

      # find which is empty - remove
      tokeep <- !(sapply(pppObjs, spatstat::is.empty))
      pppObjs <- pppObjs[tokeep]

      # concatinate metadata
      mtd <- lapply(pppObjs, function(i){i$metaData})
      mtd <- do.call("rbind", mtd)

      # superimpose
      unifiedppp <- do.call(spatstat::superimpose.ppp,
                            args = c(unname(pppObjs), list(W = spWin, check = check)))

      # attach metadata
      unifiedppp$metaData <- mtd
      class(unifiedppp) <- c(class(unifiedppp), "analytePointPattern")

      return(unifiedppp)

}
