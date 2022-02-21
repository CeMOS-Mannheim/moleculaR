#' Subset analytePointPatterns
#'
#' subset a given point patterns of type `ppp` and `analytePointPatterns` based on its `metaData`.
#'
#' @param obj  an S3 object of type `ppp` and `analytePointPatterns`.
#' @param expr logical expression indicating which points are to be kept. The expression must involve
#' the column names of the `metaData` slot). Subsetting based on `marks` is not yet supported.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` after subsetting.
#'
#' @export
#' @include manualSpatstatImport.R
#'

subsetAnalytes <- function(obj, expr){

      e <- substitute(expr)
      mtd <- obj$metaData

      # keep these idx's in the metaData table
      idxToKeep <- mtd$idx[which(eval(e, mtd, parent.frame()))]

      # subset metaData table
      mtd <- mtd[which(mtd$idx %in% idxToKeep), ]

      # call spatstat's subsetting function
      obj <- subset.ppp(obj, idx %in% idxToKeep)

      # attach metadata
      obj$metaData <- mtd
      class(obj) <- c(class(obj), "analytePointPattern")

      return(obj)

}
