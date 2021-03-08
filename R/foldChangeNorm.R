#' Normalize Intensity Based on Median Fold Change to the Mean Spectrum. 
#'
#' For every spectrum in each dataset finds median fold-change to the mean spectrum and normalizes based on that. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param meanMethod: see `?MALDIquant::mergeMassPeaks`.
#' @return
#' Filtered list of `MassPeaks` objects. `noiseLevel`is added to the `metaData` slot which 
#' encodes the median noise level for the corresponding spectrum 
#'
#' @export
#'
foldChangeNorm              = function(x, meanMethod = "mean") {
       
       
       
       #// normalization based on  ----
       # Veselkov, Kirill, et al. "BASIS: High-performance bioinformatics platform for processing of large-scale mass 
       # spectrometry imaging data in chemically augmented histology." Scientific reports 8.1 (2018): 1-11.
       
       # first find the mean mass Peaks object
       dataMean = MALDIquant::mergeMassPeaks(x, method = meanMethod)
                                          
       
       
       # for every spectrum in each tissue section find median fold-change to the mean spectrum and normalize
       
       x                    = lapply(x, FUN = function(spect) {
              
              mSpect        = MALDIquant::mass(spect)
              mMean         = MALDIquant::mass(dataMean)
              iSpect        = MALDIquant::intensity(spect)
              iMean         = MALDIquant::intensity(dataMean)
              
              idxInMean     = which(mMean %in% mSpect)
              idxInSpect    = which(mSpect %in% mMean[idxInMean])
              
              fc            = iSpect[idxInSpect] / iMean[idxInMean] 
              
              medfc         = median(fc)
              
              iSpect        = iSpect / medfc
              
              MALDIquant::intensity(spect) = iSpect
              
              # also compute noise level and add it to metadata slot
              n             = median(iSpect / MALDIquant::snr(spect))
              
              MALDIquant::metaData(spect)$medianNoise = n
              
              spect
              
       })
                                                    
       return(x)                                          

       
       
       
}

