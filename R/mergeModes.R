#' Merge positive and negative ion modes
#'
#' This function covers only the negative ion mode. 
#'
#' @param negative: 	       dataset, a list of `MassPeaks` objects for negative ion mode. 
#' @param positive: 	       dataset, a list of `MassPeaks` objects for positive ion mode. 
#' @param regionsNeg:       the regions object of the negative ion mode - not used. 
#' @param regionsPos:       the regions object of the positive ion mode - not used. 
#' @param regionName:       not used. 
#' @param numCores:         number of cores to be used within \code{parallel::mclapply}.
#' @return
#' Merged list of `MassPeaks` objects with \code{negative} being the backbone.
#'
#' @export
#'
mergeModes  = function(negative, positive, regionsNeg = NULL, regionsPos = NULL, 
                       regionName = NULL, numCores = 10, plot = FALSE) {
       
       
       
       #// let's consider the negative image as reference, we can loop through each pixel in it to find the corresponding
       #   pixel based on x-y coordinates in the positive image. 
       
       #// To make things easier, any spectrum that is present in the negative ion mode and not in the positive ion mode, 
       #   the spectrum will be not be removed. It will not have positive ion mode data but that should not be a problem. 
       
       #mergedPeakList       = negative
       #mergedRegions        = regionsNeg
    
       posCoords            = MALDIquant::coordinates(positive)
       posCoords            = apply(posCoords, 1, paste0, collapse = ",")
       
       
       negCoords            = MALDIquant::coordinates(negative)
       negCoords            = apply(negCoords, 1, paste0, collapse = ",")
       
       idxInPos             = match(x = negCoords, table = posCoords)
       
       if(FALSE) # not really important
       {
              plot(regionsNeg$polyPoints[[regionName]], type = "l", asp = 1, 
                   ylim = rev(range(regionsNeg$polyPoints[[regionName]]$y)))
              lines(regionsPos$polyPoints[[regionName]], type = "l", asp = 1, 
                    ylim = rev(range(regionsPos$polyPoints[[regionName]]$y)), col = "red")
              points(regionsNeg$allPoints[[regionName]], pch = 19, cex = 0.5)
              points(regionsPos$allPoints[[regionName]][idxInPos,], pch = "*", col = "red")

              
       }
       
       #remove pixels of the negative mode that are not in the postive mode
       # keepInNegIdx         = which(!is.na(idxInPos))
       # 
       # mergedPeakList       = negative[keepInNegIdx]
       # mergedRegions$allPoints[[regionName]] = regionsNeg$allPoints[[regionName]][keepInNegIdx]
       # mergedRegions$pointsIndx[[regionName]]= regionsNeg$pointsIndx[[regionName]][keepInNegIdx]
       # 
       
       
       
       mergedPeakList       = parallel::mclapply(X = seq_along(negative), 
                                                 mc.cores = numCores,
                                                 FUN = function(ispect) {
                                                            
                                                        if(is.na(idxInPos[ispect]))
                                                        {


                                                               #mergedRegions$allPoints[[regionName]] = mergedRegions$allPoints[[regionName]][-ispect, ]
                                                               #mergedRegions$pointsIndx[[regionName]]= mergedRegions$pointsIndx[[regionName]][-ispect]

                                                               #return(NULL)
                                                               
                                                               
                                                               MALDIquant::metaData(negative[[ispect]])$mode = rep(-1L, length(negative[[ispect]]))
                                                               
                                                              
                                                               
                                                               
                                                               return(negative[[ispect]])
                                                               

                                                        }
                                                        
                                                        
                                                        dfNeg  = data.frame(mass = negative[[ispect]]@mass,
                                                                            intensity = negative[[ispect]]@intensity,
                                                                            snr = negative[[ispect]]@snr,
                                                                            fwhm = negative[[ispect]]@metaData$fwhm,
                                                                            mode = -1L)
                                                        
                                                        dfPos  = data.frame(mass = positive[[idxInPos[ispect]]]@mass,
                                                                            intensity = positive[[idxInPos[ispect]]]@intensity,
                                                                            snr = positive[[idxInPos[ispect]]]@snr,
                                                                            fwhm = positive[[idxInPos[ispect]]]@metaData$fwhm, 
                                                                            mode = 1L)
                                                        
                                                        
                                                        df     = rbind(dfNeg, dfPos)
                                                        
                                                        
                                                        df     = df[order(df$mass), ]
                                                        
                                                        md     = negative[[ispect]]@metaData
                                                        md$fwhm= df$fwhm
                                                        md$mode= df$mode
                                                        
                                                        
                                                        MALDIquant::createMassPeaks(mass = df$mass, intensity = df$intensity, 
                                                                                    snr = df$snr, metaData = md)
                                                        
                                                            
                                                     })
       
      
       # emptyObjects         = sapply(mergedPeakList, is.null)
       # 
       # if(any(emptyObjects)){
       #        
       #        #warning("empty object(s) found in the merged dataset.")
       #        mergedPeakList= mergedPeakList[-which(emptyObjects)]
       #        
       # }

       
       # return(list(mergedPeakList = mergedPeakList,
       #             mergedRegions  = mergedRegions))
       
       return(mergedPeakList)
       
}

