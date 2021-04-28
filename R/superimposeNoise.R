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
#' @param returnMat:    if \code{TRUE} returns a sparse matrix \code{Matrix::dgCMatrix}. Otherwise, returns a list of
#' \code{MALDIquant::MassPeaks} objects (default).
#' @param verbose:      whether to print progress.
#' @return
#' A list of `MassPeaks` objects contaminated with noise, if \code{returnMat = FALSE)} or a sparse matrix of type \code{Matrix::dgCMatrix}
#' otherwise.
#'
#' @export
#'
#'

superimposeNoise        = function(x = NULL, spmat = NULL, method, mz, fwhm, noiseFactor = 1L,
                                   searchFactor = 3L, sigmaInterfering = 2L, mzTrim = 0,
                                   numSpikedPeaks = 10L, returnMat = FALSE, verbose = FALSE) {


      #// checks
      if(is.null(x) & is.null(spmat)){
            stop("one of x or spmat has to be provided. \n")
      }

      if(is.null(spmat)){

            coords      = MALDIquant::coordinates(x) # save a copy of the coordinates


            if(!MALDIquant::isMassPeaksList(x)){
                  stop("x is not a list of MassPeaks objects. \n")
            }
      }

      if(is.null(x)){

            coords      = attr(spmat, "coordinates") # save a copy of the coordinates

            if(is.null(coords))
                  stop("spmat does not hold pixel coordianetes in its attr!\n")


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
            if(verbose) {cat("creating sparse matrix representation .. \n")}
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
            fwhm  = fwhm[!is.na(idx)]
            mz    = mz[!is.na(idx)]
      }



      # add noise to data
      if(verbose) {cat("applying noise .. \n")}

      switch (method,
              "Gaussian" = {

                    #n = vector("list", length(idx))
                    n = lapply(idx, function(icol){

                        colData = .extractCol(spmat, icol)
                        #colData = spmat[ , icol]
                        abs(rnorm(nrow(spmat),
                                 mean(colData),
                                 sd(colData) * noiseFactor))

                    })

                    nmat = Matrix::sparseMatrix(i = rep(seq_len(nrow(spmat)), times = length(idx)),
                                                j = rep(idx, each = nrow(spmat)),
                                                x = unlist(n, use.names = FALSE),
                                                dimnames = list(NULL, NULL),
                                                dims = dim(spmat))


                    spmat = spmat + nmat

                    # #// to show progress
                    # pb = utils::txtProgressBar(min = 1, max = length(idx), width = 20, style = 3)
                    # cnt= 1
                    #
                    # for(i in idx){
                    #
                    #       utils::setTxtProgressBar(pb, cnt)
                    #       spmat[ , i] = spmat[ , i] + n[[i]]
                    #       cnt  = cnt + 1
                    #
                    # }
                    #
                    # close(pb)

              },
              "spiked" = {


                    idxSpike    = sample(seq(1, nrow(spmat)), numSpikedPeaks)

                    n = lapply(idx, function(icol){

                          rpois(length(idxSpike), max(.extractCol(spmat, icol), na.rm = TRUE) * 3)

                    })

                    nmat = Matrix::sparseMatrix(i = rep(idxSpike, times = length(idx)),
                                                j = rep(idx, each = length(idxSpike)),
                                                x = unlist(n, use.names = FALSE),
                                                dimnames = list(NULL, NULL),
                                                dims = dim(spmat))


                    spmat = spmat + nmat

                    # n = vector("list", ncol(spmat))
                    #
                    # n[idx] = lapply(idx, function(icol){
                    #
                    #       rpois(length(idxSpike), max(.extractCol(spmat, icol), na.rm = TRUE) * 3)
                    #
                    # })
                    #
                    # #// to show progress
                    # pb = utils::txtProgressBar(min = 1, max = length(idx), width = 20, style = 3)
                    # cnt= 1
                    #
                    # for(i in idx){
                    #
                    #       utils::setTxtProgressBar(pb, cnt)
                    #       spmat[idxSpike , i] = spmat[idxSpike , i] + n[[i]]
                    #       cnt  = cnt + 1
                    #
                    # }
                    #
                    # close(pb)

              },
              "interfering" = {

                    #n = vector("list", ncol(spmat))

                    # n[idx] = lapply(idx, function(icol){
                    #       colData = .extractCol(spmat, icol)
                    #       abs(rnorm(nrow(spmat),
                    #                 mean(colData),
                    #                 sd(colData) * noiseFactor))
                    #
                    # })

                    n = lapply(idx, function(icol){
                          colData = .extractCol(spmat, icol)
                          abs(rnorm(nrow(spmat),
                                    mean(colData),
                                    sd(colData) * noiseFactor))

                    })


                    #// to show progress
                    if(verbose & (length(idx) > 1)){
                          pb = utils::txtProgressBar(min = 1, max = length(idx), width = 20, style = 3)
                          cnt= 1
                    }


                    for(i in 1:length(idx)){


                          if(verbose & (length(idx) > 1)) {utils::setTxtProgressBar(pb, cnt); cnt  = cnt + 1}

                          # add a new mass
                          m0      = mz[i]
                          m1      = m0 + (fwhm[i] / 2.355 * sigmaInterfering)
                          idxi    = findInterval(x = m1, mzAxis)
                          mzAxis  = append(x = mzAxis, values = m1, after = idxi)

                          # m0      = mzAxis[idx[i]]
                          # m1      =  m0 + (fwhm[i] / 2.355 * sigmaInterfering)
                          # idxi    = findInterval(x = m1, mzAxis)
                          # mzAxis  = append(x = mzAxis, values = m1, after = idxi)

                          spmat   = cbind(spmat[ , seq(1,idxi)],
                                          Matrix::Matrix(data = n[[i]], nrow = nrow(spmat), sparse = TRUE),
                                          spmat[ , seq((idxi+1),ncol(spmat))])

                          colnames(spmat) = as.character(mzAxis)


                    }

                    if(verbose & (length(idx) > 1)) {close(pb)}

              },
              stop("given method is incorrect. \n"))


      if(returnMat){
            attr(spmat, "coordinates") = coords
            if(verbose) {cat("Done.\n")}
            return(spmat)
      }

      #// convert the sparse matrix to a list of MassPeaks objects
      if(verbose) {cat("re-creating MS Data as a list of MassPeaks objects .. \n")}
      x        =.spmat2MassPeaks(spmat, coords)

      if(verbose) {cat("Done.\n")}
      return(x)





}


.spmat2MassPeaks = function(spmat, coords) {

      mzAxis   = as.numeric(colnames(spmat))
      spmat    = Matrix::t(spmat) # to keep working with column-major matrix

      x        = lapply(seq_len(ncol(spmat)), function(i){

            #s     = spmat[i,] # extremely slow!
            s     = .extractCol(spmat, i)
            nnz   = which(s > 0, useNames = FALSE)

            MALDIquant::createMassPeaks(mass = mzAxis[nnz],
                                        intensity = s[nnz],
                                        metaData = list(imaging = list(pos = coords[i, ])))

      })

      return(x)
}

#// custom accessor to extract row-values from Matrix::dgCMatrix
# see https://stackoverflow.com/questions/47997184/extraction-speed-in-matrix-package-is-very-slow-compared-to-regular-matrix-class
.extractCol  = function(m, i) {
      # m: matrix of class Matrix::dgCMatrix
      # i: row index

      r       = numeric(nrow(m))   # set up zero vector for results

      # handles empty columns
      inds = seq(from=m@p[i]+1,
                 to=m@p[i+1],
                 length.out=max(0, m@p[i+1] - m@p[i]))


      r[m@i[inds]+1] = m@x[inds]     # set values

      return(r)
}
