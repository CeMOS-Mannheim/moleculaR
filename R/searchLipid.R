#' Find lipids in MSI datasets
#'
#' This function check if a deprotonated/protonated or any alkali version lipid is detectable in the MS data represented 
#' by a sparse matrix (negative or positive mode, accordingly). 
#'
#' @param m: 	       the mass to be queried (Da).
#' @param fwhm:      the fwhm at `m`.
#' @param massAxis:  the mass axis of the whole MS dataset to be searched against. 
#' @param spData:    the sparse matrix containing the negative or poisitive ion mode info. 
#' @param coords:    the coordinates of spectra contained within the spData.
#' @param metaspace: a data frame of the metaspace identifications (positive or negative mode accordingly).  
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
searchLipid          = function(m, fwhm, massAxis, spData, coords, metaspace = NA,
                                adduct = NA, mode = NA, modeAdduct = NA, 
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
       sgma          = fwhm / 2.355 # sigma i.e. std
       fiveSgma      = sgma * 5
       
       #// check if listed in metaspace
       if(!identical(metaspace, NA)) { # if metaspace identification were supplied
            
               # pool the confirmed/identified masses of metaspace 
               mtspMass      = sort(as.numeric(metaspace$mz))
               
               # match closest
               idxmtsp        = MALDIquant::match.closest(m, mtspMass, (sgma * 3), NA_integer_)
               
               # record 
               mtspConf      = ifelse(identical(idxmtsp, NA_integer_), FALSE, TRUE)
               
                  
       } else { # if metaspace identifications were not supplied flag it to FALSE
               
               mtspConf      = FALSE
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




# #// de-protonated ----
# lipTmp        = searchList$swissList[[lipidSpecies]]$`Exact m/z of [M-H]-`[i]
# 
# if(!is.na(lipTmp)) { # for example there is no de-protonated version of PC!
#        
#        
#        
#       
#        
#        # find sigma of the gaussian --> corresponds to search window
#        s             = fwhmFun(lipTmp) / 2.355 # sigma i.e. std
#        fiveSgma         = s * 5
#        
#        idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, fiveSgma, NA)
#        
#        
#        if(!is.na(idx))
#        {
#               
#               idxlwr        = MALDIquant::match.closest((lipTmp - fiveSgma), .uniqueMass, fiveSgma, idx)
#               idxupr        = MALDIquant::match.closest((lipTmp + fiveSgma), .uniqueMass, fiveSgma, idx)
#               
#               
#               coi           = as(spDataNeg[ , (idxlwr:idxupr), drop = FALSE], "matrix")  #columns of interest
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
#                                                           mode = "negative",
#                                                           modeAdduct = "M-H",
#                                                           lipidID = searchList$swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                           sumformula = searchList$swissList[[lipidSpecies]]$`Formula (pH7.3)`[i],
#                                                           fullName = searchList$swissList[[lipidSpecies]]$Name[i],
#                                                           abbrev = searchList$swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                           numDoubleBonds = searchList$swissList[[lipidSpecies]]$numDoubleBond[i],
#                                                           stringsAsFactors = F))
#               }
#               
#        }
#        
# }
# 
