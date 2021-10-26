#' Creates a sparse matrix from a list of 'MassPeaks' objects
#'
#' The sparse matrix is created via 'Matrix::sparseMatrix'. This type
#' of representation speeds-up computation compared to 'MALDIquant::MassPeaks'.
#'
#' @param x:    a list of \code{MassPeaks objects}.
#' @return
#' An S3 class 'sparseIntensityMatrix' with three slots; 'spmat' - sparse matrix
#' of type 'Matrix::dgCMatrix' holding the spectral intensities (see \code{?Matrix::sparseMatrix}),
#' 'mzAxis' - a numeric holding all m/z bins (in Da) and, finally, 'coordinates' - a two-column matrix
#' holding the coordinates of the spectra which are stored in rows of 'spmat'.
#'
#' @export
#'

createSparseMat             <- function(x) {

        if(!MALDIquant::isMassPeaksList(x)) {
                stop("Error in CreateSparseMat; x must be a list of MassPeaks objects.\n")
        }

        .mass                <- unlist(lapply(x, MALDIquant::mass))
        .intensity           <- unlist(lapply(x, MALDIquant::intensity))
        .uniqueMass          <- sort.int(unique(.mass))
        n                    <- lengths(x)
        r                    <- rep.int(seq_along(x), n)
        i                    <- findInterval(.mass, .uniqueMass)


        spmat                <- Matrix::sparseMatrix(i = r, j = i, x = .intensity,
                                                    dimnames = list(NULL, NULL),
                                                    dims = c(length(x), length(.uniqueMass)))


        return(sparseIntensityMatrix(spmat = spmat,
                                     mzAxis = .uniqueMass,
                                     coordinates = as.data.frame(MALDIquant::coordinates(x))))


}

