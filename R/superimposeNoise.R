#' Superimpose noise over MSI data
#'
#' This contaminates MSI data with additive noise at certain m/z locations. This is used
#' for testing purposes.
#'
#' @param x: 	dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param spmat: a corresponding
#' @param method: character string, the method used to add noise c("Gaussian","spiked","interfering")
#' @param mz:     numeric, m/z value where the noise will be added.
#' @param fwhm:   numeric, the fwhm at the specified m/z.
#' @param noiseFactor:    integer, this will be multiplied by the standard deviation of intensities at \code{mz} which will
#' be plugged in in \code{rnorm} to generate the intensities of the Gaussian noise.
#' @param searchFactor:    integer, search tolerance, how many standard deviations (sigma) around \code{mz} to include in the search of
#' \code{mz} peak presence.
#' @param sigmaInterfering: integer, how many standard deviations (sigma) away from \code{mz} to place the interfering peak.
#' @param mzTrim:   numeric, this will trim the m/z axis to only include \code{[mz-mzTrim,mz+mzTrim]}. Skipped if set to \code{0}
#' (default) or \code{length(mz) > 1}.
#' @param numSpikedPeaks:     integer, number of spiked peaks for \code{spiked} method.
#' @return
#' A list of `MassPeaks` objects contaminated with noise.
#'
#' @export
#'
#'

superimposeNoise        = function(x = NULL, spmat = NULL, method, mz, fwhm, noiseFactor = 3L,
                                   searchFactor = 3L, sigmaInterfering = 2L, mzTrim = 0,
                                   numSpikedPeaks = 10L) {


      #// checks
      if(is.null(x) & is.null(spmat)){
         stop("one of x or spmat has to be provided. \n")
      }

      if(is.null(spmat)){
         if(!MALDIquant::isMassPeaksList(x)){
            stop("x is not a list of MassPeaks objects. \n")
         }
      }

      if(is.null(x)){
         if(class(spmat) != "dgCMatrix"){
            stop("spmat must be of sparse matrix type Matrix::dgCMatrix. \n")
         }
      }

      if((method == "interfering") & (sigmaInterfering > searchFactor)) {
            stop("sigma interfering must be smaller than searchFactor. \n")
      }


      if((mzTrim != 0) & (length(mz) == 1) & (!is.null(spmat))){
            x       = MALDIquant::trim(x, range = c(mz-mzTrim, mz+mzTrim))
      }

      # create sparse matrix
      if(is.null(spmat)){
         spmat       = moleculaR::createSparseMat(x)
      }
      mzAxis      = as.numeric(colnames(spmat))

      # generate noise vector based on intensities of the query mass
      idx         = MALDIquant::match.closest(x = mz, table = mzAxis,
                                              tolerance = (fwhm / 2.355) * searchFactor)



      # check - mz peak has not been detected
      if(all(is.na(idx))){
            stop("No macthes for mz were found in the dataset. \n")
      }

      if(any(is.na(idx))){
            warning(paste0("m/z ", idx[is.na(idx)], " was not found in the dataset.\n"))
            idx   = idx[!is.na(idx)]
      }

      #intVec      = spmat[ , idx, drop = FALSE]


      # add noise to data
      switch (method,
             "Gaussian" = {

                      n = vector("list", ncol(spmat))
                      # n[idx] = lapply(seq_len(ncol(intVec)), function(icol){
                      n[idx] = lapply(idx, function(icol){
                            abs(rnorm(nrow(spmat),
                                      mean(spmat[ , icol]),
                                      sd(spmat[ , icol]) * noiseFactor))

                      })

                      for(i in idx){

                         spmat[ , i] = spmat[ , i] + n[[i]]

                      }


             },
             "spiked" = {

                   idxSpike    = sample(seq(1, nrow(spmat)), numSpikedPeaks)

                   n = vector("list", ncol(spmat))

                   n[idx] = lapply(idx, function(icol){

                      rpois(length(idxSpike), max(spmat[ , icol], na.rm = TRUE) * 3)

                   })


                   for(i in idx){

                      spmat[idxSpike , i] = spmat[idxSpike , i] + n[[i]]

                   }


             },
             "interfering" = {

                n = vector("list", ncol(spmat))

                n[idx] = lapply(idx, function(icol){
                   abs(rnorm(nrow(spmat),
                             mean(spmat[ , icol]),
                             sd(spmat[ , icol]) * noiseFactor))

                })

                for(i in idx){

                   # to fix: make mz as reference instead of idx

                   # add a new mass
                   m0      = mzAxis[i]
                   m1      =  m0 + (fwhm / 2.355) * sigmaInterfering
                   idxi    = findInterval(x = m1, mzAxis)
                   mzAxis  = append(x = mzAxis, values = m1, after = idxi)

                   spmat   = cbind(spmat[ , seq(1:idxi)],
                                   Matrix::Matrix(data = n[[i]], nrow = nrow(spmat), sparse = TRUE),
                                   spmat[ , seq((idxi+1):ncol(spmat))])

                   colnames(spmat) = as.character(mzAxis)

                }


             },
             stop("given method is incorrect. \n"))


      #// convert the sparse matrix to a list of MassPeaks objects
      x        =.spmat2MassPeaks(spmat)

      return(x)





}


.spmat2MassPeaks = function(spmat) {

   x        = lapply(seq_len(nrow(spmat)), function(i){

      s     = unname(spmat[i, ])
      nnz   = which(s > 0, useNames = FALSE)

      MALDIquant::createMassPeaks(mass = as.numeric(colnames(spmat)[nnz]),
                                  intensity = s[nnz],
                                  metaData = list(imaging = list(pos = attr(spmat, "coordinates")[i, ])))

   })

   return(x)
}

