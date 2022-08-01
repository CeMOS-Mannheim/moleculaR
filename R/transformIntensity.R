#' Intensity transformation for SPPs
#'
#'
#'
#' @param spp: 	  a spatial point pattern of type 'ppp' with its `marks` being a dataframe with
#' two columns `idx`and `intensity`.
#' @param method:  a character specifying the method of transformation, one of
#' `c("z-score", "sqrt", "scaling")` such that `z-score` represents (analyte-wise)
#' standardization (subtraction of mean and division by std), `sqrt` is the square
#' root transformation and `scaling` is linear mapping of intensities into [0,1] range.
#' @param irange:    a numeric vector of two elements speciying the output range of
#' the scalled intensities.
#' @param forceRange: logical, whether to force-map the intensities in to the range
#' given by `irange` regardless of the `method` used.
#'
#' @return
#'
#'
#' @export
#' @include manualSpatstatImport.R
#'

transformIntensity <- function(spp, method = "z-score", irange = c(0, 1), forceRange = TRUE){

         if(!("analytePointPattern" %in% class(spp))){
                  stop("spp must be of 'analytePointPattern' and 'ppp' class. \n")
         }

         if(spp$n == 0){
                  warning("Provided 'spp' has zero points, transformation not possible.\n")
                  return(spp)
         }

         if(diff(irange) <= 0){
                  stop("incorrect 'irange' provided. \n")
         }

          # offset the lowest value (zero) because it is reserved to pixels with empty points in spatstat
         if(min(irange) == 0){	
              irange[1] <- irange[1] + 0.0001
         }

         spp <- switch(method,

                       "z-score" = .zscore(spp),
                       "sqrt" = .sqrt(spp),
                       "scaling" = .scale(spp, irange),
                       stop("Unrecognized 'method', must be one of",
                            " c('z-score', 'sqrt', 'scaling').\n")
         )

         if(forceRange & method != "scaling"){
                  spp <- .scale(spp, irange)
         }

         return(spp)

}

.zscore <- function(spp){

         numAnalytes <- length(spp$metaData$mzVals)

         for(ii in 1 : length(numAnalytes)){

                  # which idx
                  idx <- which(spp$marks$idx == spp$metaData$idx[ii])


                  spp$marks$intensity[idx] <- (spp$marks$intensity[idx] - mean(spp$marks$intensity[idx])) /
                           sd(spp$marks$intensity[idx])

         }

         return(spp)

}

.sqrt <- function(spp){

         spp$marks$intensity <- sqrt(spp$marks$intensity)

         return(spp)

}

.scale <- function(spp, irange){

         numAnalytes <- length(spp$metaData$mzVals)

         for(ii in 1 : length(numAnalytes)){

                  # which idx
                  idx <- which(spp$marks$idx == spp$metaData$idx[ii])

                  spp$marks$intensity[idx] <- .linMap(spp$marks$intensity[idx], irange[1], irange[2])

         }

         return(spp)
}

#.linMap <- function(x, irange) approxfun(range(x), irange)(x)
# linear mapping of values into a custom range
.linMap <- function(i, a, b) approxfun(range(i), c(a, b))(i)
