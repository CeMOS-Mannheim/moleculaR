#' Superimpose noise over MSI data
#'
#' This contaminates MSI data with additive noise at certain m/z locations. This is used
#' for testing purposes.
#'
#' @param x: 	dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param spData: a corresponding sparse MSI matrix of type 'sparseIntensityMatrix'.
#' @param method: character string, the method used to add noise c("Poisson", "Gaussian","spiked","interfering").
#' @param mz:     numeric, m/z value where the noise will be added.
#' @param fwhm:   numeric, the fwhm at the specified m/z.
#' @param noiseFactor:    integer, this will be multiplied by the standard deviation of intensities at \code{mz} which will
#' be plugged into \code{rnorm} to generate the intensities of the Gaussian noise.
#' @param searchFactor:    integer, search tolerance, how many standard deviations (sigma) around \code{mz} to include in the search of
#' \code{mz} peak presence.
#' @param sigmaInterfering: integer, how many standard deviations (sigma) away from \code{mz} to place the interfering peak.
#' @param mzTrim:   numeric, this will trim the m/z axis to only include \code{[mz-mzTrim,mz+mzTrim]}. Skipped if set to \code{0}
#' (default) or \code{length(mz) > 1}.
#' @param numSpikedPeaks:     integer, number of spiked peaks for \code{spiked} method.
#' @param returnMat:    if \code{TRUE} returns a sparse MSI matrix of type `sparseIntensityMatrix`. Otherwise, returns a list of
#' \code{MALDIquant::MassPeaks} objects (default).
#' @param verbose:      whether to print progress.
#' @return
#' A list of `MassPeaks` objects contaminated with noise, if \code{returnMat = FALSE)} or a sparse matrix of type \code{Matrix::dgCMatrix}
#' otherwise.
#'
#' @export
#'
#'

superimposeNoise        = function(x = NULL, spData = NULL, method, mz, fwhm, noiseFactor = 1L,
                                   searchFactor = 3L, sigmaInterfering = 2L, mzTrim = 0,
                                   numSpikedPeaks = 10L, returnMat = FALSE, verbose = FALSE) {


         #// checks
         if(is.null(x) & is.null(spData)){
                  stop("one of x or spData has to be provided. \n")
         }

         if(is.null(spData)){

                  if(!MALDIquant::isMassPeaksList(x)){
                           stop("x is not a list of MassPeaks objects. \n")
                  }

                  coords      = MALDIquant::coordinates(x) # save a copy of the coordinates

         }

         if(is.null(x)){

                  coords      = spData$coordinates # save a copy of the coordinates

                  if(is.null(coords))
                           stop("spData does not hold pixel coordianetes in its attr!\n")


                  if(class(spData) != "sparseIntensityMatrix"){
                           stop("spData must be of type 'sparseIntensityMatrix'. \n")
                  }
         }

         if((method == "interfering") & (sigmaInterfering > searchFactor)) {
                  stop("sigma interfering must be smaller than searchFactor. \n")
         }


         if((mzTrim != 0) & (length(mz) == 1) & (!is.null(spData))){
                  x       = MALDIquant::trim(x, range = c(mz-mzTrim, mz+mzTrim))
         }

         # create sparse matrix
         if(is.null(spData)){
                  if(verbose) {cat("creating sparse matrix representation .. \n")}
                  spData       = createSparseMat(x)
         }

         mzAxis      = spData$mzAxis

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
                 "Poisson" = {

                          n = lapply(idx, function(icol){

                                   colData = .extractCol(spData$spmat, icol)
                                   abs(rpois(nrow(spData$spmat),
                                             mean(colData)) * noiseFactor)

                          })

                          nmat = Matrix::sparseMatrix(i = rep(seq_len(nrow(spData$spmat)), times = length(idx)),
                                                      j = rep(idx, each = nrow(spData$spmat)),
                                                      x = unlist(n, use.names = FALSE),
                                                      dimnames = list(NULL, NULL),
                                                      dims = dim(spData$spmat))


                          spData$spmat = spData$spmat + nmat


                 },
                 "Gaussian" = {

                          n = lapply(idx, function(icol){

                                   colData = .extractCol(spData$spmat, icol)
                                   abs(rnorm(nrow(spData$spmat),
                                             mean(colData),
                                             sd(colData)) * noiseFactor)

                          })

                          nmat = Matrix::sparseMatrix(i = rep(seq_len(nrow(spData$spmat)), times = length(idx)),
                                                      j = rep(idx, each = nrow(spData$spmat)),
                                                      x = unlist(n, use.names = FALSE),
                                                      dimnames = list(NULL, NULL),
                                                      dims = dim(spData$spmat))


                          spData$spmat = spData$spmat + nmat


                 },
                 "spiked" = {


                          idxSpike    = sample(seq(1, nrow(spData$spmat)), numSpikedPeaks)

                          n = lapply(idx, function(icol){

                                   rpois(length(idxSpike), max(.extractCol(spData$spmat, icol), na.rm = TRUE) * 3)

                          })

                          nmat = Matrix::sparseMatrix(i = rep(idxSpike, times = length(idx)),
                                                      j = rep(idx, each = length(idxSpike)),
                                                      x = unlist(n, use.names = FALSE),
                                                      dimnames = list(NULL, NULL),
                                                      dims = dim(spData$spmat))


                          spData$spmat = spData$spmat + nmat



                 },
                 "interfering" = {



                          n = lapply(idx, function(icol){

                                   colData = .extractCol(spData$spmat, icol)
                                   abs(rpois(nrow(spData$spmat),
                                             mean(colData)) * noiseFactor)

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


                                   spData$spmat   = cbind(spData$spmat[ , seq(1,idxi)],
                                                          Matrix::Matrix(data = n[[i]], nrow = nrow(spData$spmat), sparse = TRUE),
                                                          spData$spmat[ , seq((idxi+1),ncol(spData$spmat))])

                                   spData$mzAxis = mzAxis

                          }

                          if(verbose & (length(idx) > 1)) {close(pb)}

                 },
                 stop("given method is incorrect. \n"))


         if(returnMat){

                  if(verbose) {cat("Done.\n")}

                  return(spData)
         }

         #// convert the sparse matrix to a list of MassPeaks objects
         if(verbose) {cat("re-creating MS Data as a list of MassPeaks objects .. \n")}
         x        = .spmat2MassPeaks(spData)

         if(verbose) {cat("Done.\n")}

         return(x)





}


.spmat2MassPeaks = function(spData) {

         tmpMat    = Matrix::t(spData$spmat) # to keep working with column-major matrix

         x        = lapply(seq_len(ncol(tmpMat)), function(i){

                  #s     = spmat[i,] # extremely slow!
                  s     = .extractCol(tmpMat, i)
                  nnz   = which(s > 0, useNames = FALSE)

                  MALDIquant::createMassPeaks(mass = spData$mzAxis[nnz],
                                              intensity = s[nnz],
                                              metaData = list(imaging = list(pos = spData$coordinates[i, ])))

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
