#' Exports datasets to imzML for metaspace validation
#'
#' The export is done via `MALDIquantForeign`. The resolution at the highest peak is appended to `fileName`. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param fileName: A custom name for the file, without extention. 
#' @param path: path where to write the imzML file. 
#' @param pixelSize: Yes, you guessed it, Pixel Size!
#' @return
#' Nothing. check the path for the newly created files.   
#'
#' @export
#'
metaSpaceExport      = function(x, fileName = "MSIdata", path, pixelSize = c(50, 50)) {
       
       
       
       #// export the regions to metaspace for verification ----
       
       # compute resolution at a specific peak
       mxSnr         = lapply(x, FUN = function(spect) {
              
              mxIdx  = which.max(spect@snr)
              
              data.frame(mass = spect@mass[mxIdx], peakRes = spect@mass[mxIdx] / spect@metaData$fwhm[mxIdx])
              
       })
       
       mxSnr         = do.call("rbind", mxSnr)
       mxSnr$mass    = round(mxSnr$mass, 6)
       mostFreqMass  = names(sort(table(mxSnr$mass), decreasing = TRUE)[1])
       massRes       = mean(mxSnr$peakRes[mxSnr$mass == as.numeric(mostFreqMass)])
       
       
       # convert to massspectrum objects
       x             = parallel::mclapply(x, mc.cores = 1, FUN = function(spect) {
              
                     MALDIquant::createMassSpectrum(mass = spect@mass, intensity = spect@intensity, metaData = spect@metaData)
              
                     })
       
       #export as imzml
       MALDIquantForeign::exportImzMl(x = x, 
                                      path = file.path(path, paste0(fileName,"-res ", round(massRes), " at mass ", mostFreqMass, ".imzML")),
                                      processed = TRUE, 
                                      coordinates = MALDIquant::coordinates(x),
                                      pixelSize = pixelSize)
       
       return(NULL)
       
       
}

