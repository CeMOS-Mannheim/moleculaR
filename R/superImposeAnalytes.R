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

superimposeAnalytes <- function(pppObjs, spWin = NULL, check = TRUE){

      # find which is empty - remove
      # tokeep <- !(sapply(pppObjs, function(i) {
      #       spatstat.geom::is.empty.ppp(i) | is.null(i)
      #       }))

      # torm <- sapply(pppObjs, function(i) {
      #       (spatstat.geom::is.empty.ppp(i) | is.null(i))
      # })
      #
      # if(all(torm)){
      #       if(is.null(spWin)){
      #             return(spatstat.geom::ppp(x = integer(0), y = integer(0)))
      #       } else{
      #             return(spatstat.geom::ppp(x = integer(0), y = integer(0), window = spWin))
      #       }
      # }
      #
      # pppObjs <- pppObjs[-torm]

      # concatinate metadata
      mtd <- lapply(pppObjs, function(i){i$metaData})
      mtd <- do.call("rbind", mtd)

      # superimpose
      unifiedppp <- do.call(spatstat.geom::superimpose.ppp,
                            args = c(unname(pppObjs), list(W = spWin, check = check)))

      # attach metadata
      unifiedppp$metaData <- mtd
      class(unifiedppp) <- c(class(unifiedppp), "analytePointPattern")

      return(unifiedppp)

}
