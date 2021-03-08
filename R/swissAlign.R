#' Blind alignment of peaks to theoretical masses in the Swisslipids DB
#'
#' This function covers only the negative ion mode. 
#'
#' @param x: 	Dataset, a list of `MassPeaks` objects. 
#' @param swissdf: The Swisslipids database loaded as a dataframe. 
#' @param tol: relative tolerance.  
#' @return
#' Re-aligned peaks list. 
#'
#' @export
#'
swissAlign  = function(x, swissdb, tol = 2e-6) {
       
       
       
       # #// try re-aligning the MassPeaks objects with respect to all lipids found in swissdblipids ----

       #   first create the reference MassPeaks object
       # Note: PC/PE/SM/LPC/DAG -> positive ion mode and PI/PE/PS/PA/PG -> negative ion mode.
       
       piIdx                = grep(pattern = "PI", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       peIdx                = grep(pattern = "PE", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5882520/
       psIdx                = grep(pattern = "PS", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       paIdx                = grep(pattern = "PA", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       pgIdx                = grep(pattern = "PG", x = swissdb$`Abbreviation*`, ignore.case = FALSE)

       swissdbNeg             = swissdb[c(piIdx, peIdx, psIdx, paIdx, pgIdx), ]


       refPeak              = MALDIquant::createMassPeaks(mass = sort(swissdbNeg$`Exact m/z of [M-H]-`),
                                                          intensity = rep(1L, nrow(swissdbNeg)),
                                                          snr = rep(3L, nrow(swissdbNeg)))


       wf                   = MALDIquant::determineWarpingFunctions(l = x, reference = refPeak,
                                                 tolerance = tol, method="lowess",
                                                 allowNoMatches=TRUE, plot=FALSE,
                                                 plotInteractive = FALSE)





       x                    = MALDIquant::warpMassPeaks(l = x, w = wf, emptyNoMatches = FALSE)


       
       
       return(x)
       
}

