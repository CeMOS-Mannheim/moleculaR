#' Load msi data 
#'
#' This is to load msi data stored in peaks.sqlite DBs. For this two files are required; .mis and .sqlite. 
#' The \code{regName} defines which region (within the measurement run) is to be loaded and processed. 
#'
#' @param sqlitePath:       path to the sqlite file containing centroided data. 
#' @param misPath:          path to the mis file.
#' @param FunLibPath:       path to the supporting functions library of this project. 
#' @param regName:          a string indicating the name of the region to be loaded and processed (as appears in the mis file). 
#' @return
#' A list of msData and the correspoinding region. 
#' 
#' @export
#'
loadData                    = function(sqlitePath, misPath,  regName, ionMode) {
       
       
       
       #// load the dataset ----
       cat("loading FTICR dataset .. \n")
       peaksList             = msiTools::GetCentroidedFTICR(sqlitePath = sqlitePath,
                                                           misPath = misPath,
                                                           plotRegions = FALSE, includeManualRegions = TRUE)
       
       # two element list; $peaksList and $regions
     
       
       
       #// find tissue regions ----
       regIdxNeg            = grep(x=peaksList$regions$minMaxCoords$ROI_name, 
                                   pattern = paste0("(", regName, ")", collapse = "|"))   #"[0-9]"
       
       #// focus on one region
       peaksList$peaksList  = peaksList$peaksList[unlist(peaksList$regions$pointsIndx[regIdxNeg])]
       

       
       #// compute the fwhm and the corresponding relative tolerance to use it for binning
       cat("computing fwhm .. \n")

       fwhmFun              = estimateFwhm(x = peaksList$peaksList, plot = FALSE, savePlot = FALSE)
       
       
       #// bin peaks  ----
       cat("peak binning .. \n")
       peaksList$peaksList  =  MALDIquant::binPeaks(peaksList$peaksList, 
                                                       tolerance = fwhmFun(400)/400, #focusing on lipids 
                                                       method = "relaxed")
                                                    
       

       
       
       #// filter out peaks which occur in less than 1% of the time - the built-in function of MALDIquant produces errors ----
       cat("peaks filtering .. \n")
       peaksList$peaksList = filterPeaks(x = peaksList$peaksList, minFreq = 0.01, ionMode = ionMode)
              
              

       #// normalization based on  ----
       # Veselkov, Kirill, et al. "BASIS: High-performance bioinformatics platform for processing of large-scale mass 
       # spectrometry imaging data in chemically augmented histology." Scientific reports 8.1 (2018): 1-11.
       cat("intensity normalization .. \n")
       peaksList$peaksList = foldChangeNorm(peaksList$peaksList)
                                                
       cat("Done. \n")
       

       return(list(msData = peaksList$peaksList, 
                   regions = peaksList$regions, 
                   regName = regName, 
                   coordinates = MALDIquant::coordinates(peaksList$peaksList),
                   fwhmFun = fwhmFun))
       
       
       
}


