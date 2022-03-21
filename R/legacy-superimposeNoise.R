#' Superimpose noise over MSI data - only for testing
#'
#' This contaminates MSI data with additive noise at certain m/z locations. This is used
#' for testing purposes.
#'
#' @param x: 	dataset, a list of \code{MALDIquant::MassPeaks} objects.
#' @param method: character string, the method used to add noise c("Gaussian","spiked","interfering")
#' @param mz:     numeric, m/z value where the noise will be added.
#' @param fwhm:   numeric, the fwhm at the specified m/z.
#' @param noiseFactor:    integer, this will be multiplied by the standard deviation of intensities at \code{mz} which will
#' be plugged in in \code{rnorm} to generate the intensities of the Gaussian noise.
#' @param searchFactor:    integer, search tolerance, how many standard deviations (sigma) around \code{mz} to include in the search of
#' \code{mz} peak presence.
#' @param sigmaInterfering: integer, how many standard deviations (sigma) away from \code{mz} to place the interfering peak.
#' @param mzTrim:   numeric, this will trim the m/z axis to only include \code{[mz-mzTrim,mz+mzTrim]}.
#' Defaults to 10 and if set to 0 or \code{length(mz) > 1}, trimming will be skipped.
#' @param numSpikedPeaks:     integer, number of spiked peaks for \code{spiked} method.
#' @return
#' A list of `MassPeaks` objects contaminated with noise.
#'
#'
#'
#'

superimposeNoise        = function(x, method, mz, fwhm, noiseFactor = 3L,
                                   searchFactor = 3L, sigmaInterfering = 2L, mzTrim = 10L,
                                   numSpikedPeaks = 10L) {

      #// checks
      if(!MALDIquant::isMassPeaksList(x)){
            stop("x is not a list of MassPeaks objects. \n")
      }

      if((method == "interfering") & (sigmaInterfering > searchFactor)) {
            stop("sigma interfering must be smaller than searchFactor. \n")
      }

      if((mzTrim != 0) & (length(mz) == 1)){
          x       = MALDIquant::trim(x, range = c(mz-mzTrim, mz+mzTrim))
      }


      # create sparse matrix
      sptmp       = moleculaR::createSparseMat(x)

      # generate noise vector based on intensities of the query mass
      idx         = MALDIquant::match.closest(x = mz, table = as.numeric(colnames(sptmp)),
                                              tolerance = (fwhm / (2*sqrt(2*log(2)))) * searchFactor)



      # check - mz peak has not been detected
      if(all(is.na(idx))){
            stop("No macthes for mz were found in the dataset. \n")
      }

      if(any(is.na(idx))){
            warning(paste0("m/z ", idx[is.na(idx)], " was not found in the dataset.\n"))
            idx   = idx[!is.na(idx)]
      }

      intVec      = sptmp[ , idx]


      # add noise to data
      x           = switch (method,
                            "Gaussian" = {
                                  n = rnorm(length(intVec), mean(intVec), sd(intVec) * noiseFactor)

                                  mapply(FUN = .gaussian,
                                         x = x, n = n,
                                         SIMPLIFY = FALSE, USE.NAMES = FALSE,
                                         MoreArgs = list(idx = idx))
                            },
                            "spiked" = {
                                  idxSpike    = sample(seq(1, length(x)), numSpikedPeaks)
                                  n           = rep(0, length(x))
                                  n[idxSpike] = rpois(length(idxSpike), max(intVec, na.rm = TRUE)*3)

                                  mapply(FUN = .gaussian,
                                         x = x, n = n,
                                         SIMPLIFY = FALSE, USE.NAMES = FALSE,
                                         MoreArgs = list(idx = idx))
                            },
                            "interfering" = {
                                  n = rnorm(length(intVec), mean(intVec), sd(intVec) * noiseFactor)

                                  mapply(FUN = .interfering,
                                         x = x, n = n,
                                         SIMPLIFY = FALSE, USE.NAMES = FALSE,
                                         MoreArgs = list(mz = mz, fwhm = fwhm,
                                                         meanNoise = mean(intVec),
                                                         sdNoise = sd(intVec) * noiseFactor,
                                                         sigmaInterfering = sigmaInterfering))
                            },
                            stop("given method is incorrect. \n"))

      return(x)

}

.gaussian   = function(x, n, idx) {

      x@intensity[idx] = x@intensity[idx] + n

      x


}



.interfering   = function(x, n, mz, fwhm, sigmaInterfering) {

      m = x@mass
      s = x@intensity
      snr= x@snr
      md = x@metaData

      # add a new mass
      m1 = mz + (fwhm / (2*sqrt(2*log(2)))) * sigmaInterfering

      m  = append(x = m, values = m1, after = idx)
      s = append(x = s, values = n, after = idx)
      snr = append(x = snr, values = mean(snr), after = idx)

      MALDIquant::createMassPeaks(mass = m, intensity = s, snr = snr, metaData = md)


}



# superimposeNoise        = function(x, method, mz, fwhm, noiseFactor = 3L,
#                                    searchFactor = 3L, sigmaInterfering = 2L, mzTrim = 10L,
#                                    numSpikedPeaks = 10L) {
#
#    #// checks
#    if(!MALDIquant::isMassPeaksList(x)){
#       stop("x is not a list of MassPeaks objects. \n")
#    }
#
#    if((method == "interfering") & (sigmaInterfering > searchFactor)) {
#       stop("sigma interfering must be smaller than searchFactor. \n")
#    }
#
#    if((mzTrim != 0) & (length(mz) == 1)){
#       x       = MALDIquant::trim(x, range = c(mz-mzTrim, mz+mzTrim))
#    }
#
#
#    # create sparse matrix
#    sptmp       = moleculaR::createSparseMat(x)
#
#    # generate noise vector based on intensities of the query mass
#    idx         = MALDIquant::match.closest(x = mz, table = as.numeric(colnames(sptmp)),
#                                            tolerance = (fwhm / 2.355) * searchFactor)
#
#
#
#    # check - mz peak has not been detected
#    if(all(is.na(idx))){
#       stop("No macthes for mz were found in the dataset. \n")
#    }
#
#    if(any(is.na(idx))){
#       warning(paste0("m/z ", idx[is.na(idx)], " was not found in the dataset.\n"))
#       idx   = idx[!is.na(idx)]
#    }
#
#    intVec      = sptmp[ , idx, drop = FALSE]
#
#
#    # add noise to data
#    x           = switch (method,
#                          "Gaussian" = {
#
#                             n = lapply(seq_len(length(x)), function(i){ # noise list whose length = length(x)
#                                sapply(seq_len(ncol(intVec), function(icol){ # each list entry should have a vector of
#                                   # noise values at the given mz.
#
#                                   abs(rnorm(1,  # noise value cannot be negative
#                                             mean(intVec[ ,icol]),
#                                             sd(intVec[ ,icol]) * noiseFactor))
#
#                                })
#                             })
#
#                                mapply(FUN = .gaussian,
#                                       xi = x, ni = n,
#                                       SIMPLIFY = FALSE, USE.NAMES = FALSE,
#                                       MoreArgs = list(idx = idx))
#
#                          },
#                          "spiked" = {
#
#                             idxSpike    = sample(seq(1, length(x)), numSpikedPeaks)
#
#                             n = lapply(seq_len(length(x)), function(i){ # noise list whose length = length(x)
#
#                                if(!(i %in% idxSpike)){ # if i not in idxSpike return a zero-valued vector of n
#                                   return(rep(0, length(idx)))
#                                }
#
#                                sapply(seq_len(ncol(intVec), function(icol){ # each list entry should have a vector of
#                                   # noise values at the given mz.
#
#                                   rpois(1, max(intVec[ , icol], na.rm = TRUE) * 3)
#
#
#                                })
#                             })
#
#
#                                mapply(FUN = .gaussian,
#                                       x = x, n = n,
#                                       SIMPLIFY = FALSE, USE.NAMES = FALSE,
#                                       MoreArgs = list(idx = idx))
#                          },
#                          "interfering" = {
#
#                             n = lapply(seq_len(length(x)), function(i){ # noise list whose length = length(x)
#                                sapply(seq_len(ncol(intVec), function(icol){ # each list entry should have a vector of
#                                   # noise values at the given mz.
#
#                                   abs(rnorm(1,  # noise value cannot be negative
#                                             mean(intVec[ ,icol]),
#                                             sd(intVec[ ,icol]) * noiseFactor))
#
#                                })
#                             })
#
#                                mapply(FUN = .interfering,
#                                       x = x, n = n,
#                                       SIMPLIFY = FALSE, USE.NAMES = FALSE,
#                                       MoreArgs = list(mz = mz, fwhm = fwhm, idx = idx,
#                                                       meanNoise = mean(intVec),
#                                                       sdNoise = sd(intVec) * noiseFactor,
#                                                       sigmaInterfering = sigmaInterfering))
#
#                          },
#                          stop("given method is incorrect. \n"))
#
#                             return(x)
#
# }
#
# .gaussian   = function(xi, ni, idx) {
#
#
#
#    xi@intensity[idx] = xi@intensity[idx] + ni
#
#    xi
#
#
# }
#
#
#
# .interfering   = function(xi, ni, mz, fwhm, sigmaInterfering) {
#
#    m = xi@mass
#    s = xi@intensity
#    snr= xi@snr
#    md = xi@metaData
#
#    # add a new mass
#    m1 = mz + (fwhm / 2.355) * sigmaInterfering
#
#    m  = append(x = m, values = m1, after = idx)
#    s = append(x = s, values = ni, after = idx) ## multi mz case?
#    snr = append(x = snr, values = mean(snr), after = idx)
#
#    MALDIquant::createMassPeaks(mass = m, intensity = s, snr = snr, metaData = md)
#
#
# }
