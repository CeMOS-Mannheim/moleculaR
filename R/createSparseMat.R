#' Creates a sparse matrix from a list of `MassPeaks` objects
#'
#' The sparse matrix is created via `Matrix::sparseMatrix`. Two attributes are added to the resulting sparse matrix; `featuresMode` which 
#' encodes the polarity (negative or postive ion mode) in which a peak was detected (if both modes are combined into one dataset) and `noiseLevel`. 
#' Note that in order to have these two attributes you have to run `filterPeaks` and `foldChangeNorm` beforehand. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @return
#' A sparse matrix, see `?Matrix::sparseMatrix` for more info. 
#'
#' @export
#'
createSparseMat             = function(x) {
       
       
       
       
       # first create a spars matrix per tissue region
       .mass                = unlist(lapply(x, MALDIquant::mass))
       .mode                = unlist(lapply(x, FUN = function(spect) spect@metaData$mode))
       .intensity           = unlist(lapply(x, MALDIquant::intensity))
       .noise               = unlist(lapply(x, FUN = function(spect) spect@metaData$medianNoise))
       
       if(is.null(.mode)) {
              .mode         = rep(1L, length(.mass))
              warning("ionMode is absent from the metadata and is ignored. See FilterPeaks to include it.")
       }
       
       if(is.null(.noise)) {
              .noise        = rep(0, length(x))
              warning("noiseLevel is absent from the metadata and is ignored. See foldChangeNorm to include it.")
       }
       
       
       .signedMass          = .mass * .mode #to later know the "mode" origin of masses
       .unsortedUniqueMass  = unique(.signedMass)
       .tmp                 = sort.int(abs(.unsortedUniqueMass), index.return = TRUE)
       
       .uniqueMass          = .tmp$x
       .uniqueMode          = sign(.unsortedUniqueMass[.tmp$ix])
       
       n                    = lengths(x)
       r                    = rep.int(seq_along(x), n)
       i                    = findInterval(.mass, .uniqueMass)
       
       spmat                = Matrix::sparseMatrix(i = r, j = i, x = .intensity, 
                                                   dimnames = list(NULL, .uniqueMass),
                                                   dims = c(length(x), length(.uniqueMass)))
       
       attr(spmat, "featuresMode") = .uniqueMode
       attr(spmat, "noiseLevel") = .noise
       
       
       
       return(spmat)                                          
       
       
       
       
}

