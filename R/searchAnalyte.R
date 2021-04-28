#' Find an analyte in the dataset
#'
#' This function check if an analyte (given by m/z) is detectable in the MS data represented 
#' by a sparse matrix
#'
#' @param m: 	the mass to be queried (Da).
#' @param fwhm: the fwhm at `m`.
#' @param massAxis: the mass axis of the whole MS dataset to be searched against. 
#' @param spData: the sparse matrix containing the data. 
#' @param coords: the coordinates of spectra contained within the spData.
#' @param wMethod: wighting method; c("Gaussian", "sum", "max", "mean"). 
#' @return
#' A dataframe containing the hits for `m` in `msData`. 
#'
#' @export
#'
searchAnalyte     = function(m, fwhm, massAxis, spData, coords, wMethod = "Gaussian") {
       
       
       df            = data.frame(x = integer(0), 
                                  y = integer(0),
                                  mass = numeric(0),
                                  intensity = numeric(0), 
                                  stringsAsFactors = F)
       
       

       # find sigma of the gaussian --> corresponds to search window
       s             = fwhm / 2.355 # sigma i.e. std
       fiveS         = s * 5
       
       idx           = MALDIquant::match.closest(m , massAxis, fiveS, NA)
       
       
       if(!is.na(idx))
       {
              
              idxlwr        = MALDIquant::match.closest((m - fiveS), massAxis, fiveS, idx)
              idxupr        = MALDIquant::match.closest((m + fiveS), massAxis, fiveS, idx)
              
              
              coi           = as(spData[ , (idxlwr:idxupr), drop = FALSE], "matrix")  #columns of interest

              combinedCols  = switch(wMethod, 
                     "Gaussian" = {
                            gw = gaussWeight(x = massAxis[(idxlwr:idxupr)], 
                                          m = m, 
                                          fwhm = fwhm)
              
                            coi           = sweep(coi, MARGIN = 2, gw, "*")
              
                            rowSums(coi)

                     },
                     "sum" = {
                            rowSums(coi)
                     },
                     "mean" = {
                            rowMeans(coi)
                     },
                     "max" = {
                            apply(coi, 1, max)
                     })
              
              
              
              
              detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
              
              if(length(detectedIn) > 0) {
                     
                     # coordinates of these spectra
                     detectedCoord = coords[detectedIn, , drop = F]
                     
                     df            = rbind(df, data.frame(x = detectedCoord[ , "x"],
                                                          y = detectedCoord[ , "y"],
                                                          mass = m,
                                                          intensity = combinedCols[detectedIn],
                                                          
                                                          stringsAsFactors = F))
              }
              
       }
       
       
       
       
       return(df)                                          
       
       
       
       
}



