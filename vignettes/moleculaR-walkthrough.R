## ----conf, include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----setup, eval=FALSE, include=TRUE------------------------------------------
#  # not run
#  install.packages("devtools")
#  devtools::install_github("CeMOS-Mannheim/moleculaR", build_vignettes=TRUE)
#  
#  library(moleculaR)
#  

## ----importData, eval=FALSE, include=TRUE-------------------------------------
#  #-- not run --#
#  
#  # read MSI data
#  imzmlFile         <- "pathToFile.imzML"
#  msData            <- readCentrData(path = imzmlFile)
#  
#  # single spectrum
#  spectrFile        <- "pathToFile.tsv"
#  msSpectr          <- readSingleSpect(spectrFile)
#  

## ----import-mtspc, eval=FALSE, include=TRUE-----------------------------------
#  #-- not run --#
#  
#  # load the metaspace annotations file
#  pathToMtspc       <- "pathToFile.csv"
#  mtspc             <- read.csv(file = pathToMtspc, skip = 2, header = TRUE, colClasses = "character")
#  

## ----import-sldb, eval=FALSE, include=TRUE------------------------------------
#  #-- not run --#
#  
#  # load the processed swisslipids db
#  pathTosldb        <- system.file("extdata", "swisslipids-speciesOnly-sep2020.tsv", package = "moleculaR", mustWork = TRUE)
#  sldb              <- loadSwissDB(pathTosldb)
#  

## ----fwhm-esitmation, eval=FALSE, include=TRUE--------------------------------
#  #-- not run --#
#  
#  # estimate fwhm from msSpectr
#  fwhmObj           <- estimateFwhm(s = msSpectr)
#  

## ----fwhm-defaults, eval=FALSE----------------------------------------------
#  # Approximate instrument models when no profile spectrum is available
#  fwhmObj <- estimateFwhm() # FTICR: resolving power 123000 at m/z 249
#  fwhmObj <- estimateFwhm(modality = "TOF-reflector")
#  fwhmObj <- estimateFwhm(modality = "TIMS-TOF")
#  fwhmObj <- estimateFwhm(userMassRes = c(massRes = 150000, mass = 400))
#  getFwhm(fwhmObj, 400) # 400 / 150000
#  # Measured widths remain preferable: estimateFwhm(s = msSpectr)

## ----preproc, eval=FALSE, include=TRUE----------------------------------------
#  #-- not run --#
#  
#  # bin peaks
#  msData            <-  MALDIquant::binPeaks(msData,
#                                             tolerance = getFwhm(fwhmObj, 400)/400, #focusing on lipids
#                                             method = "relaxed")
#  
#  
#  # filter out peaks which occur in less than 1% of the time
#  # Note: use 'moleculaR::' namespace to distinguish it from MALDIquant::filterPeaks if MALDIquant is loaded.
#  msData            <- moleculaR::filterPeaks(x = msData, minFreq = 0.01)
#  
#  

## ----load---------------------------------------------------------------------
library(moleculaR)

data("processed-example-Data")

# to see the loaded objects
ls()


## ----plotFwhm-----------------------------------------------------------------
plot(fwhmObj)

## ----getFwhm------------------------------------------------------------------
# FWHM at m/z 400
getFwhm(fwhmObj, 400)

## ----sparseMat----------------------------------------------------------------
#// create sparse matrix representation
spData <- createSparseMat(x = msData)

## ----query--------------------------------------------------------------------
# input by m/z value
queryMass         <- 788.5447 

## ----ionImage-----------------------------------------------------------------
# compute the regular ion image - returns an AnalytePiontPattern
sppIonImage      <- searchAnalyte(m = queryMass, 
                                  fwhm = getFwhm(fwhmObj, queryMass), 
                                  spData = spData, 
                                  wMethod = "sum")


# compute a raster image of the sppIonImage 
ionImage        <- spatstat.geom::pixellate(sppIonImage,
                          weights = sppIonImage$marks$intensity,
                          W = spatstat.geom::as.mask(sppIonImage$window,
                                                dimyx=c(diff(sppIonImage$window$yrange),
                                                        diff(sppIonImage$window$xrange))),
                          padzero = FALSE)



## ----sppmoi, results='hide', warning=FALSE------------------------------------
# compute spatial point pattern of the analyte
sppMoi          <- searchAnalyte(m = queryMass, 
                            fwhm = getFwhm(fwhmObj, queryMass), 
                            spData = spData,
                            wMethod = "Gaussian")


plotAnalyte(sppMoi, analyte = paste0("m/z ", round(queryMass, 4)))

## ----MPM, results='hide', fig.width=12, fig.height=16, warning=FALSE----------
#// compute MPM - default parameters
probImg         <- probMap(sppMoi)

txt  <- paste0("m/z ", round(queryMass, 4), " ± ", round(getFwhm(fwhmObj, queryMass), 4))

#// plot everything together
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)
plot(probImg, what = "detailed", analyte = txt, ionImage = ionImage)



## ----ion-image-matrix, eval=FALSE-------------------------------------------
#  myMatrix <- matrix(NA_real_, nrow = 20, ncol = 15)
#  myMatrix[3:18, 3:13] <- 1
#  myMatrix[7:13, 6:10] <- 10
#  sppFromMatrix <- ionImage2spp(myMatrix, orientation = "yx")
#  mpmFromMatrix <- probMap(sppFromMatrix, bw = 1, seed = 1)
#  plot(mpmFromMatrix, what = "MPM")
#  # Optional mass metadata:
#  sppFromMatrix <- ionImage2spp(myMatrix, orientation = "yx", mzVal = 788.5447)

## ----ion-image-mpm, eval=FALSE----------------------------------------------
#  # MALDIquant is an additional input format, using [x, y, mass].
#  ionArray <- MALDIquant::msiSlices(msData, center = queryMass,
#                                    tolerance = 0.05)
#  sppFromImage <- ionImage2spp(ionArray)
#  mpmFromImage <- probMap(sppFromImage, seed = 1)
#  plot(mpmFromImage, what = "MPM")
#  ionArrays <- MALDIquant::msiSlices(msData, center = c(queryMass, 810.6),
#                                     tolerance = 0.05)
#  sppSecond <- ionImage2spp(ionArrays, slice = 2)

## ----CPPMs, results='hide'----------------------------------------------------
cat("Batch lipid search is ongoing - this will take several minutes - \n")

lipidHits <- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb,
                               adduct = c("M+H", "M+Na", "M+K"),
                               numCores = 4, 
                               verifiedMasses = as.numeric(mtspc$mz),
                               confirmedOnly = TRUE, verbose = TRUE)

lipidHits

## ----metadata-----------------------------------------------------------------
# show metaData 
head(lipidHits$metaData)

## ----lipidClasses-------------------------------------------------------------
# whow all detected lipid classes
unique(lipidHits$metaData$lipidClass)

## ----CPPM.plot, results='hide', fig.width=12, fig.height=16, warning=FALSE----
# choose one lipid species
lipidClass <- "PA(x:x)"

# subset lipidHits
paHits <- subsetAnalytes(lipidHits, lipidClass == "PA(x:x)")


# compute MPM - default parameters
probImg    <- probMap(paHits, bwMethod = "scott", sqrtTansform = TRUE)


#// plot everything together
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)
plot(probImg, what = "detailed", 
      analyte = paste0(lipidClass, " - n=", 
                       length(probImg$sppMoi$metaData$mzVals)))


## ----lyso-GPLs----------------------------------------------------------------
# create a subset representing lysoGPLs
ofInterest <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)","LPI(x:x)", "LPS(x:x)", 
                "PA(x:x)", "PC(x:x)", "PE(x:x)","PG(x:x)", "PI(x:x)", "PS(x:x)")

# subset lipidHits
lysoGPLs <- subsetAnalytes(lipidHits, lipidClass %in% ofInterest)

lysoGPLs

# detected classes
cat("detected lipid classes: \n")
unique(lysoGPLs$metaData$lipidClass)

## ----metadataCols-------------------------------------------------------------
# meta data table columns
colnames(lysoGPLs$metaData)

## ----adducts,fig.width=12, fig.height=10--------------------------------------
# detected adducts
cat("detected adducts: ")
table(lysoGPLs$metaData$adduct)


# subsetting
kHits           <- subsetAnalytes(lysoGPLs, adduct == "M+K")
# or NHits      <- subsetAnalytes(lysoGPLs, adduct == "M+Na")
# or HHits      <- subsetAnalytes(lysoGPLs, adduct == "M+H")


## ----ions.plot.K, results='hide', fig.width=12, fig.height=16, warning=FALSE----

# compute CPPM
probImg    = probMap(kHits, bwMethod = "scott", sqrtTansform = TRUE)

if(probImg$sppMoi$n > 50000) {
  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
}
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)

plot(probImg, what = "detailed", analyte = paste0("[M+k]+ of (lyso)GPLs - n=", length(probImg$sppMoi$metaData$mzVals)))




## ----subsetting.saturation----------------------------------------------------
detectedSaturation <- c("sat", "mono-unsat", "di-unsat", "poly-unsat")


# subsetting
satHits             <- subsetAnalytes(lysoGPLs, numDoubleBonds == 0)
# or monoHits       <- subsetAnalytes(lysoGPLs, numDoubleBonds == 1)
# or diHits         <- subsetAnalytes(lysoGPLs, numDoubleBonds == 2)
# or polyHits       <- subsetAnalytes(lysoGPLs, numDoubleBonds > 2)


## ----cppms.unsat.plot,results='hide', fig.width=12, fig.height=16, warning=FALSE----
# compute CPPM
probImg    = probMap(satHits, bwMethod = "scott", sqrtTansform = TRUE)

if(probImg$sppMoi$n > 50000) {
  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
}
par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)
plot(probImg, what = "detailed", analyte = paste0("mono-unsat of ", "(lyso)GPLs", " - n=", length(probImg$sppMoi$metaData$mzVals)))


## -----------------------------------------------------------------------------
sessionInfo()

