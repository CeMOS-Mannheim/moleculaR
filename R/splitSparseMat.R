#' Split a sparse matrix into negative and positive ion modes
#'
#' This function splits a sparse matrix containing the spectral information into 
#' negative and positive ion mode sparse matrices. 
#'
#' @param x: 	The original sparse matrix containing both negative and positive ion modes. 
#' It has to contain the attribute `featureMode`. 
#' @return
#' A list of two sparse matrices, one for each mode. 
#'
#' @export
#'
splitSparseMat             = function(x) {
       
       
       
       
       # the table of masses to be searched against
       #.uniqueMass          = as.numeric(x@Dimnames[[2]])
       .ionMode             = attr(x, "featuresMode")
       if(is.null(.ionMode)){
              stop("the featuresMode attribute of the sparse matrix x is empty.\n")
       }
       
       spdims               = dim(x)
       
       # create two sparse matrecis; one for each mode
       idxNeg               = which(.ionMode < 0)
       spDataNeg            = x
       spDataNegMask        = Matrix::sparseMatrix(i = rep(seq(1, spdims[1]), each = length(idxNeg)),
                                                   j = rep(idxNeg, spdims[1]),
                                                   x = 1L,
                                                   dims = spdims)
       
       spDataNeg            = spDataNeg * spDataNegMask
       #invisible(gc())
       
       
       idxPos               = which(.ionMode > 0)
       spDataPos            = x
       spDataPosMask        = Matrix::sparseMatrix(i = rep(seq(1, spdims[1]), each = length(idxPos)),
                                                   j = rep(idxPos, spdims[1]),
                                                   x = 1L,
                                                   dims = spdims)
       
       spDataPos            = spDataPos * spDataPosMask
       invisible(gc())
       
       
       
       return(list(spDataNeg = spDataNeg, 
                   spDataPos = spDataPos))                                          
       
       
       
       
}

