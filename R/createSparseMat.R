#' Creates a sparse matrix from a list of \code{MassPeaks} objects
#'
#' The sparse matrix is created via \code{Matrix::sparseMatrix}.
#'
#' @param x:    a list of \code{MassPeaks objects}.
#' @return
#' A sparse matrix, see \code{?Matrix::sparseMatrix} for more info.
#'
#' @export
#' @keywords internal
#'
createSparseMat             = function(x) {

        if(!MALDIquant::isMassPeaksList(x)) {
                stop("Error in CreateSparseMat; x must be a list of MassPeaks objects.\n")
        }

        .mass                = unlist(lapply(x, MALDIquant::mass))
        .intensity           = unlist(lapply(x, MALDIquant::intensity))
        .uniqueMass          = sort.int(unique(.mass))
        n                    = lengths(x)
        r                    = rep.int(seq_along(x), n)
        i                    = findInterval(.mass, .uniqueMass)

        # spmat                = Matrix::sparseMatrix(i = r, j = i, x = .intensity,
        #                                             dimnames = list(NULL, .uniqueMass),
        #                                             dims = c(length(x), length(.uniqueMass)))
        #
        spmat                = Matrix::sparseMatrix(i = r, j = i, x = .intensity,
                                                    dimnames = list(NULL, NULL),
                                                    dims = c(length(x), length(.uniqueMass)))



       return(new("sparseIntensityMatrix",
                  i = spmat@i,
                  p = spmat@p,
                  Dim = spmat@Dim,
                  Dimnames = spmat@Dimnames,
                  x = spmat@x,
                  factors = spmat@factors,
                  mzAxis = .uniqueMass,
                  coordinates = MALDIquant::coordinates(x)))

}

