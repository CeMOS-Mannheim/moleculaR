#' Filter Peaks
#'
#' This function filters out peaks which occur in less than `minFreq` percent of all pixels. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param minFreq: Peaks occuring in `<= minFreq * length(x)` will be omitted. 
#' @param ionMode: Either `+1L` or `-1L` indicating the polarity with which the dataset was measured. 
#' @return
#' Filtered list of `MassPeaks` objects. `mode` is added to the `metaData` slot;  which 
#' encodes the polarity (negative or postive ion mode) in which a peak was detected (useful 
#' if both modes are combined into one dataset).
#'
#' @export
#'
filterPeaks  = function(x, minFreq = 0.01, ionMode = -1L) {
       
       
       
       #// filter out peaks which occur in less than 1% of the time - the built-in function of MALDIquant produces errors ----
       
       # first create a spars matrix per tissue region
       .mass                = unlist(lapply(x, MALDIquant::mass))
       .intensity           = unlist(lapply(x, MALDIquant::intensity))
       .uniqueMass          = sort.int(unique(.mass))
       n                    = lengths(x)
       r                    = rep.int(seq_along(x), n)
       i                    = findInterval(.mass, .uniqueMass)
       
       spmat                = Matrix::sparseMatrix(i = r, j = i, x = .intensity, 
                                                   dimnames = list(NULL, .uniqueMass),
                                                   dims = c(length(x), length(.uniqueMass)))
       
       
       
       # find which peaks (columns) occur in less than 1% of the tissue
       relFreq              = diff(spmat@p) / nrow(spmat)
       kpIdx                = which(relFreq >= 0.01)
       kpMass               = .uniqueMass[kpIdx]
       
       x                    = lapply(x, FUN = function(i) {
              
              xm            = MALDIquant::mass(i)
              xkpIdx        = which(xm %in% kpMass)
              mt            = MALDIquant::metaData(i)
              mt$fwhm       = mt$fwhm[xkpIdx]
              mt$mode       = rep(ionMode, length(xm[xkpIdx])) # just for compatibiltiy with previous notebooks -> shows peak ion mode
              
              MALDIquant::createMassPeaks(mass = xm[xkpIdx], 
                                          intensity = MALDIquant::intensity(i)[xkpIdx],
                                          snr = MALDIquant::snr(i)[xkpIdx],
                                          metaData = mt)
              
       })
       

       
       
       return(x)
       
}

