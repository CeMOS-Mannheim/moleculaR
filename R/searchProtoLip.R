#' Find protonated version of a lipid
#'
#' This function check if a protonated version lipid is detectable in the MS data represented 
#' by a sparse matrix (positive mode). 
#'
#' @param m: 	the mass to be queried (Da).
#' @param fwhm: the fwhm at `m`.
#' @param massAxis: the mass axis of the whole MS dataset to be searched against. 
#' @param spData: the sparse matrix containing the positive ion mode info. 
#' @param coords: the coordinates of spectra contained within the spData.
#' @param lipidID
#' @param sumformula
#' @param fullName
#' @param abbrev
#' @param numDoubleBonds
#' 
#' @return
#' A dataframe containing the hits for `m` in `msData`. 
#'
#' @export
#'
searchProtoLip     = function(m, fwhm, massAxis, spData, coords, 
                                lipidID = NA, sumformula = NA, fullName = NA,
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
                                  fullName = character(0),
                                  abbrev = character(0),
                                  numDoubleBonds = integer(0),
                                  stringsAsFactors = F)
       
       
       
       #// de-protonated ----
       
       
       
       # find sigma of the gaussian --> corresponds to search window
       s             = fwhm / 2.355 # sigma i.e. std
       fiveS         = s * 5
       
       idx           = MALDIquant::match.closest(m , massAxis, fiveS, NA)
       
       
       if(!is.na(idx))
       {
              
              idxlwr        = MALDIquant::match.closest((m - fiveS), massAxis, fiveS, idx)
              idxupr        = MALDIquant::match.closest((m + fiveS), massAxis, fiveS, idx)
              
              
              coi           = as(spData[ , (idxlwr:idxupr), drop = FALSE], "matrix")  #columns of interest
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
                                                          adduct = "H1",
                                                          mode = "positive",
                                                          modeAdduct = "M+H",
                                                          lipidID = lipidID,
                                                          sumformula = sumformula,
                                                          fullName = fullName,
                                                          abbrev = abbrev,
                                                          numDoubleBonds = numDoubleBonds,
                                                          stringsAsFactors = F))
              }
              
       }
       
       
       
       
       return(df)                                          
       
       
       
       
}




# #// protonated ----
# lipTmp        = searchList$swissList[[lipidSpecies]]$`Exact m/z of [M+H]+`[i]
# 
# if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species
#        
#        
#        # find sigma of the gaussian --> corresponds to search window
#        s             = fwhmFun(lipTmp) / 2.355 # sigma i.e. std
#        fiveS         = s * 5
#        
#        idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, fiveS, NA)
#        
#        
#        if(!is.na(idx))
#        {
#               
#               idxlwr        = MALDIquant::match.closest((lipTmp - fiveS), .uniqueMass, fiveS, idx)
#               idxupr        = MALDIquant::match.closest((lipTmp + fiveS), .uniqueMass, fiveS, idx)
#               
#               
#               coi           = as(spDataPos[ , (idxlwr:idxupr), drop = FALSE], "matrix")  #columns of interest
#               gw            = gaussWeight(x = .uniqueMass[(idxlwr:idxupr)], 
#                                           m = lipTmp, 
#                                           fwhm = fwhmFun(lipTmp))
#               
#               coi           = sweep(coi, MARGIN = 2, gw, "*")
#               
#               combinedCols  = rowSums(coi)
#               
#               
#               
#               detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#               
#               if(length(detectedIn) > 0) {
#                      
#                      # coordinates of these spectra
#                      detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]
#                      
#                      df            = rbind(df, data.frame(x = detectedCoord[ , "x"],
#                                                           y = detectedCoord[ , "y"],
#                                                           intensity = combinedCols[detectedIn],
#                                                           mass = lipTmp,
#                                                           adduct = "H1",
#                                                           mode = "positive",
#                                                           modeAdduct = "M+H",
#                                                           lipidID = searchList$swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                           sumformula = searchList$swissList[[lipidSpecies]]$`Formula (pH7.3)`[i],
#                                                           fullName = searchList$swissList[[lipidSpecies]]$Name[i],
#                                                           abbrev = searchList$swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                           numDoubleBonds = searchList$swissList[[lipidSpecies]]$numDoubleBond[i],
#                                                           stringsAsFactors = F))
#               }
#        }
#        
# }