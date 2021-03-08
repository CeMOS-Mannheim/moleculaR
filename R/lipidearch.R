# # loop through all lipids 
# for(lipidSpecies in names(searchList$lipidHits))
# {
#        
#        
#        tic                  = Sys.time()
#        
#        #lipidSpecies = names(searchList$lipidHits)[1]
#        #i = 1
#        
#        searchList$lipidHits[[lipidSpecies]] = parallel::mclapply(X = seq(1, nrow(swissList[[lipidSpecies]])), mc.cores = 1, FUN = function(i) {
#               
#               
#               
#               df            = data.frame(x = integer(0), 
#                                          y = integer(0),
#                                          mass = numeric(0),
#                                          intensity = numeric(0), 
#                                          adduct = character(0), 
#                                          mode = character(0),
#                                          modeAdduct = character(0),
#                                          lipidID = character(0),
#                                          sumformula = character(0), 
#                                          fullName = character(0),
#                                          abbrev = character(0),
#                                          stringsAsFactors = F)
#               
#               #// de-protonated ----
#               lipTmp        = searchList$swissList[[lipidSpecies]]$`Exact m/z of [M-H]-`[i]
#               
#               if(!is.na(lipTmp)) { # for example there is no de-protonated version of PC!
#                      
#                      
#                      # Check if that lipid was detected in our datasets
#                      .uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
#                      
#                      idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, lipTmp * 5e-6, NA)
#                      
#                      
#                      if(!is.na(idx))
#                      {
#                             
#                             idxlwr        = MALDIquant::match.closest((lipTmp - (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             idxupr        = MALDIquant::match.closest((lipTmp + (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             
#                             if(idxlwr != idxupr) {
#                                    
#                                    combinedCols = Matrix::rowSums(spmatList[[regName]][ , (idxlwr:idxupr)])
#                                    
#                             } else {
#                                    
#                                    combinedCols = spmatList[[regName]][ , (idx)]
#                             }
#                             
#                             detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#                             
#                             # coordinates of these spectra
#                             detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]              
#                             
#                             df            = rbind(df, data.frame(x = detectedCoord[ , "x"], 
#                                                                  y = detectedCoord[ , "y"], 
#                                                                  intensity = combinedCols[detectedIn],
#                                                                  mass = lipTmp,
#                                                                  adduct = "H1",
#                                                                  mode = "negative", 
#                                                                  modeAdduct = "M-H",
#                                                                  lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                                  sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i], 
#                                                                  fullName = swissList[[lipidSpecies]]$Name[i],
#                                                                  abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                                  stringsAsFactors = F))
#                             
#                      }
#                      
#               }
#               
#               #// protonated ----
#               lipTmp        = swissList[[lipidSpecies]]$`Exact m/z of [M+H]+`[i]
#               
#               if(!is.na(lipTmp)) {
#                      
#                      # Check if that lipid was detected in our datasets
#                      #.uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
#                      
#                      idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, lipTmp * 5e-6, NA)
#                      
#                      
#                      if(!is.na(idx))
#                      {
#                             
#                             idxlwr        = MALDIquant::match.closest((lipTmp - (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             idxupr        = MALDIquant::match.closest((lipTmp + (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             
#                             if(idxlwr != idxupr) {
#                                    
#                                    combinedCols = Matrix::rowSums(spmatList[[regName]][ , (idxlwr:idxupr)])
#                                    
#                             } else {
#                                    
#                                    combinedCols = spmatList[[regName]][ , (idx)]
#                             }
#                             
#                             detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#                             
#                             # coordinates of these spectra
#                             detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]              
#                             
#                             df            = rbind(df, data.frame(x = detectedCoord[ , "x"], 
#                                                                  y = detectedCoord[ , "y"], 
#                                                                  intensity = combinedCols[detectedIn], 
#                                                                  mass = lipTmp,
#                                                                  adduct = "H1", 
#                                                                  mode = "positive",
#                                                                  modeAdduct = "M+H",
#                                                                  lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                                  sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i], 
#                                                                  fullName = swissList[[lipidSpecies]]$Name[i],
#                                                                  abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                                  stringsAsFactors = F))
#                             
#                      }
#               }
#               
#               
#               #// Na+ adduct ----
#               lipTmp        = swissList[[lipidSpecies]]$`Exact m/z of [M+Na]+`[i]
#               
#               if(!is.na(lipTmp)) {
#                      
#                      # Check if that lipid was detected in our datasets
#                      #.uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
#                      
#                      idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, lipTmp * 5e-6, NA)
#                      
#                      
#                      if(!is.na(idx))
#                      {
#                             
#                             idxlwr        = MALDIquant::match.closest((lipTmp - (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             idxupr        = MALDIquant::match.closest((lipTmp + (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             
#                             if(idxlwr != idxupr) {
#                                    
#                                    combinedCols = Matrix::rowSums(spmatList[[regName]][ , (idxlwr:idxupr)])
#                                    
#                             } else {
#                                    
#                                    combinedCols = spmatList[[regName]][ , (idx)]
#                             }
#                             
#                             
#                             detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#                             
#                             # coordinates of these spectra
#                             detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]              
#                             
#                             df            = rbind(df, data.frame(x = detectedCoord[ , "x"], 
#                                                                  y = detectedCoord[ , "y"], 
#                                                                  intensity = combinedCols[detectedIn], 
#                                                                  mass = lipTmp,
#                                                                  adduct = "Na1", 
#                                                                  mode = "positive",
#                                                                  modeAdduct = "M+Na",
#                                                                  lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                                  sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i], 
#                                                                  fullName = swissList[[lipidSpecies]]$Name[i],
#                                                                  abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                                  stringsAsFactors = F))
#                             
#                      }
#                      
#               }
#               
#               #// K+ adduct ----
#               lipTmp        = swissList[[lipidSpecies]]$`Exact m/z of [M+K]+`[i]
#               
#               if(!is.na(lipTmp)) {
#                      
#                      # Check if that lipid was detected in our datasets
#                      #.uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
#                      
#                      idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, lipTmp * 5e-6, NA)
#                      
#                      
#                      if(!is.na(idx))
#                      {
#                             
#                             idxlwr        = MALDIquant::match.closest((lipTmp - (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             idxupr        = MALDIquant::match.closest((lipTmp + (lipTmp * 2.5e-6)), .uniqueMass, lipTmp * 2.5e-6, idx)
#                             
#                             if(idxlwr != idxupr) {
#                                    
#                                    combinedCols = Matrix::rowSums(spmatList[[regName]][ , (idxlwr:idxupr)])
#                                    
#                             } else {
#                                    
#                                    combinedCols = spmatList[[regName]][ , (idx)]
#                             }
#                             
#                             detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#                             
#                             # coordinates of these spectra
#                             detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]              
#                             
#                             df            = rbind(df, data.frame(x = detectedCoord[ , "x"], 
#                                                                  y = detectedCoord[ , "y"], 
#                                                                  intensity = combinedCols[detectedIn], 
#                                                                  mass = lipTmp,
#                                                                  adduct = "K1", 
#                                                                  mode = "positive",
#                                                                  modeAdduct = "M+K",
#                                                                  lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                                  sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i], 
#                                                                  fullName = swissList[[lipidSpecies]]$Name[i],
#                                                                  abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                                  stringsAsFactors = F))
#                             
#                      }
#               }
#               
#               return(df)
#        })
#        
#        
#        
#        
#        toc                  =  Sys.time() - tic
#        
#        cat("##### Region ", regName, " #####\n")
#        
#        print(toc)
#        
#        cat("\nNumber of", lipidSpecies,"lipids detected within the region = ", 
#            length(which(!unlist(lapply(lipidHits[[lipidSpecies]], FUN = function(x) nrow(x) == 0)))),
#            " out of ", nrow(swissList[[lipidSpecies]]), lipidSpecies ," lipids within SwissLipids database. \n")
#        
#        isDetected           = !unlist(lapply(lipidHits[[lipidSpecies]], FUN = function(x) nrow(x) == 0))
#        numDetection         = unlist(lapply(lipidHits[[lipidSpecies]], FUN = function(x) nrow(x)))
#        
#        isSignificant        = numDetection > (length(e$msDataPeakList[[regName]]) * 0.01) 
#        
#        cat("Number of lipids that can be significant [# of detections > 0.01] = ", 
#            length(which(isSignificant)), "\n\n")
#        
#        
#        
#        lipidHits[[lipidSpecies]]     = do.call("rbind", lipidHits[[lipidSpecies]])
#        
#        
#        
# }





# #// de-protonated ----
# lipTmp        = searchList$swissList[[lipidSpecies]]$`Exact m/z of [M-H]-`[i]
# 
# if(!is.na(lipTmp)) { # for example there is no de-protonated version of PC!
# 
# 
#        # Check if that lipid was detected in our datasets
#        .uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
# 
#        # compute the gaussian weights for the entire row length
#        spdims        = dim(spmatList[[regName]])
#        gwVec         = gaussWeight(x = .uniqueMass, m = lipTmp,
#                                    fwhm = fwhmFun(lipTmp),
#                                    ionMode = attr(spmatList[[regName]], "featuresMode"),
#                                    plot = FALSE)
# 
#        nzi           = which(gwVec != 0) # gives the column index of non zero indecies
# 
#        # create a sparse matrix of the same dimentions as the spmatList[[regName]] only holding
#        # the gaussian weights at the specified locations
#        gwsp          = Matrix::sparseMatrix(i = rep(seq(1, spdims[1]), each = length(nzi)),
#                                             j = rep(nzi, spdims[1]),
#                                             x = gwVec[nzi],
#                                             dims = spdims)
# 
# 
#        # switch off positive ion mode signals to focus on negative ion mode
#        idx           = which(attr(spmatList[[regName]], "featuresMode") < 0)
#        spDataNegFilter= Matrix::sparseMatrix(i = rep(seq(1, spdims[1]), each = length(nzi)),
#                                              j = rep(nzi, spdims[1]),
#                                              x = gwVec[nzi],
#                                              dims = spdims)
# 
# 
#        # now multiply the gaussian weighted sp matrix with the data matrix
#        spData        = spmatList[[regName]] * gwsp
# 
# 
#        # sum the rows
#        collapsedCols = Matrix::rowSums(spData)
# 
# 
#        # which rows (i.e. vector elements) are non-zero
#        detectedIn = which(collapsedCols > 0)
# 
#        if(length(detectedIn) > 0) # there are pixels (rows) with non-zero elements (signal detected)
#        {
# 
# 
#               # coordinates of these spectra
#               detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]
# 
#               df            = rbind(df, data.frame(x = detectedCoord[ , "x"],
#                                                    y = detectedCoord[ , "y"],
#                                                    intensity = collapsedCols[detectedIn],
#                                                    mass = lipTmp,
#                                                    adduct = "H1",
#                                                    mode = "negative",
#                                                    modeAdduct = "M-H",
#                                                    lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                                    sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i],
#                                                    fullName = swissList[[lipidSpecies]]$Name[i],
#                                                    abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                                    stringsAsFactors = F))
# 
#        }
# 
# }




# #// de-protonated ----
# lipTmp        = searchList$swissList[[lipidSpecies]]$`Exact m/z of [M-H]-`[i]
# 
# if(!is.na(lipTmp)) { # for example there is no de-protonated version of PC!
#        
#        
#        # Check if that lipid was detected in our datasets
#        .uniqueMass   = as.numeric(spmatList[[regName]]@Dimnames[[2]])
#        
#        
#        #idx           = MALDIquant::match.closest(lipTmp , .uniqueMass, (3 * s), NA)
#        
#        idxlwr        = MALDIquant::match.closest((lipTmp - (s * 5)), .uniqueMass)
#        idxupr        = MALDIquant::match.closest((lipTmp + (s * 5)), .uniqueMass)
#        
#        if(idxlwr != idxupr) {
#               
#               # subset the spmatrix
#               spData = spmatList[[regName]][ , (idxlwr:idxupr)]
#               
#               
#               
#               combinedCols = Matrix::rowSums(spmatList[[regName]][ , (idxlwr:idxupr)])
#               
#               #here create the sparse or dense matrices. 
#               
#               
#        } else {
#               
#               combinedCols = spmatList[[regName]][ , (idx)]
#        }
#        
#        detectedIn    = which(combinedCols > 0) # in which spectra it had a non-zero value
#        
#        # coordinates of these spectra
#        detectedCoord = MALDIquant::coordinates(e$msDataPeakList[[regName]])[detectedIn, , drop = F]
#        
#        df            = rbind(df, data.frame(x = detectedCoord[ , "x"],
#                                             y = detectedCoord[ , "y"],
#                                             intensity = combinedCols[detectedIn],
#                                             mass = lipTmp,
#                                             adduct = "H1",
#                                             mode = "negative",
#                                             modeAdduct = "M-H",
#                                             lipidID = swissList[[lipidSpecies]]$`Lipid ID`[i],
#                                             sumformula = swissList[[lipidSpecies]]$`Formula (pH7.3)`[i],
#                                             fullName = swissList[[lipidSpecies]]$Name[i],
#                                             abbrev = swissList[[lipidSpecies]]$`Abbreviation*`[i],
#                                             stringsAsFactors = F))
#        
#        
#        
# }