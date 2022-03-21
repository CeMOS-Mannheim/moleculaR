#' Find lipids in MSI datasets
#'
#' This function check if a deprotonated/protonated or any alkali version lipid is detectable in the MS data represented
#' by a sparse matrix.
#'
#' @param m: 	     the mass to be queried (Da).
#' @param fwhm:      the fwhm at `m`.
#' @param massAxis:  the mass axis of the whole MS dataset to be searched against.
#' @param spData:    the sparse matrix containing the negative or poisitive ion mode info.
#' @param coords:    the coordinates of spectra contained within the spData.
#' @param mtspc:     a data frame of the metaspcae identifications.
#' @param confirmedOnly: if \code{TRUE}, returns detections only if confirmed by metaspace.
#' @param lipidID:   for internal use only.
#' @param sumformula for internal use only.
#' @param fullName   for internal use only.
#' @param abbrev     for internal use only.
#' @param numDoubleBonds for internal use only.
#'
#' @return
#' A dataframe containing the hits for `m` in `msData`.
#'
#'
#'
searchLipid          = function(m, fwhm, massAxis, spData, coords, mtspc = NA,
                                confirmedOnly = FALSE, adduct = NA, mode = NA, modeAdduct = NA,
                                lipidID = NA, sumformula = NA,
                                abbrev = NA, numDoubleBonds = NA) {


       df            = data.frame(x = integer(0),
                                  y = integer(0),
                                  mass = numeric(0),
                                  intensity = numeric(0),
                                  adduct = character(0),
                                  mode = character(0),
                                  modeAdduct = character(0),
                                  lipidID = character(0),
                                  sumformula = character(0),
                                  abbrev = character(0),
                                  numDoubleBonds = integer(0),
                                  confirmed = logical(0),
                                  stringsAsFactors = F)


       #// find sigma of the gaussian --> corresponds to search window
       sgma          = fwhm / (2*sqrt(2*log(2))) # sigma i.e. std
       fiveSgma      = sgma * 5

       #// check if listed in mtspc
       if(!identical(mtspc, NA)) { # if mtspc identification were supplied

               # pool the confirmed/identified masses of mtspc
               mtspMass      = sort(as.numeric(mtspc$mz))

               # match closest
               idxmtsp        = MALDIquant::match.closest(m, mtspMass, (sgma * 3), NA_integer_)

               # record
               mtspConf      = ifelse(identical(idxmtsp, NA_integer_), FALSE, TRUE)


       } else { # if mtspc identifications were not supplied flag it to FALSE

               mtspConf      = FALSE
       }

      if(confirmedOnly){ # if only confirmed detections are needed

            if(!mtspConf){ # and there are no confirmed detections then return empty df
                  return(df)
            }
      }

       #// find mass against msi dataset
       idx           = MALDIquant::match.closest(m , massAxis, fiveSgma, NA)


       if(!is.na(idx))
       {

              idxlwr        = MALDIquant::match.closest((m - fiveSgma), massAxis, fiveSgma, idx)
              idxupr        = MALDIquant::match.closest((m + fiveSgma), massAxis, fiveSgma, idx)


              coi           = Matrix::as.matrix(spData[ , (idxlwr:idxupr), drop = FALSE])  #columns of interest
              gw            = gaussWeight(x = massAxis[(idxlwr:idxupr)],
                                          m = m,
                                          fwhm = fwhm)

              coi           = sweep(coi, MARGIN = 2, gw, "*")

              combinedCols  = rowSums(coi)



              detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value

              if(length(detectedIn) > 0) {

                     # coordinates of these spectra
                     detectedCoord = coords[detectedIn, , drop = F]

                     df            = rbind(df, data.frame(x = detectedCoord[ , "x"],
                                                          y = detectedCoord[ , "y"],
                                                          intensity = combinedCols[detectedIn],
                                                          mass = m,
                                                          adduct = adduct,
                                                          mode = mode,
                                                          modeAdduct = modeAdduct,
                                                          lipidID = lipidID,
                                                          sumformula = sumformula,
                                                          abbrev = abbrev,
                                                          numDoubleBonds = numDoubleBonds,
                                                          confirmed = mtspConf,
                                                          stringsAsFactors = F))
              }

       }




       return(df)




}




