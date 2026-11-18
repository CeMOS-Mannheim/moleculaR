library(moleculaR)

near <- function(x, y, tolerance = 1e-12) {
      isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance))
}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
messages <- function(expr) {
      result <- character()
      withCallingHandlers(force(expr), message = function(m) {
            result <<- c(result, conditionMessage(m))
            invokeRestart("muffleMessage")
      })
      result
}

# Omitting s and explicit NULL select the same calibrated FTICR default.
notice <- messages(default <- estimateFwhm())
explicit <- suppressMessages(estimateFwhm(s = NULL))
masses <- c(150, 249, 400, 800, 1100)
stopifnot(inherits(default, "fwhm"),
          near(getFwhm(default, masses), getFwhm(explicit, masses)),
          near(getFwhm(default, 249), 249 / 123000),
          length(notice) == 1L,
          grepl("FTICR", notice), grepl("123000", notice), grepl("249", notice))

for(modality in c("FTICR", "TOF-reflector", "TIMS-TOF", "TOF", "timsTOF")){
      for(anchor in list(NULL, c(30000, 600), c(mass = 400, massRes = 150000))){
            model <- suppressMessages(estimateFwhm(modality = modality, userMassRes = anchor))
            reference <- useDefaultFwhm(modality = modality, userMassRes = anchor)
            stopifnot(near(getFwhm(model, masses), getFwhm(reference, masses)))
      }
}
notice <- messages(custom <- estimateFwhm(modality = "timsTOF",
                                          userMassRes = c(mass = 400, massRes = 150000)))
stopifnot(grepl("TIMS-TOF", notice), grepl("150000", notice), grepl("400", notice))
fixed <- suppressMessages(estimateFwhm(modality = "TOF", constant = TRUE,
                                      userMassRes = c(30000, 600)))
stopifnot(near(getFwhm(fixed, masses), rep(600 / 30000, length(masses))))

# Input errors must not silently turn into approximate models.
for(bad in list(123, "spectrum", list(), NA, FALSE)){
      notice <- messages(stopifnot(fails(estimateFwhm(s = bad))))
      stopifnot(length(notice) == 0L)
}
stopifnot(fails(estimateFwhm(modality = "unknown")),
          fails(estimateFwhm(userMassRes = c(0, 249))))

# A measured profile and a list of profiles still recover Gaussian peak widths.
mass <- seq(100, 200, by = 0.002)
centres <- seq(110, 190, by = 4)
intensity <- 1 + Reduce(`+`, lapply(centres, function(m) 100 * exp(-0.5 * ((mass - m) / 0.05)^2)))
spectrum <- MALDIquant::createMassSpectrum(mass = mass, intensity = intensity)
notice <- messages(measured <- estimateFwhm(spectrum))
measuredList <- estimateFwhm(list(spectrum), spectraSampling = 1)
expected <- 2 * sqrt(2 * log(2)) * 0.05
stopifnot(length(notice) == 0L, inherits(measured, "fwhm"),
          near(getFwhm(measured, centres), rep(expected, length(centres)), tolerance = 0.03),
          near(getFwhm(measuredList, centres), getFwhm(measured, centres)))
centroided <- MALDIquant::createMassPeaks(mass = centres,
      intensity = rep(100, length(centres)), snr = rep(10, length(centres)),
      metaData = list(fwhm = rep(expected, length(centres))))
precomputed <- estimateFwhm(list(centroided), spectraSampling = 1)
stopifnot(near(getFwhm(precomputed, centres), rep(expected, length(centres))))
flat <- MALDIquant::createMassSpectrum(mass = mass, intensity = rep(1, length(mass)))
stopifnot(fails(estimateFwhm(flat)))

# Plotting and SVG export work for both measured and approximate models.
pdfPath <- tempfile(fileext = ".pdf")
svgPath <- tempfile(fileext = ".svg")
grDevices::pdf(pdfPath)
device <- grDevices::dev.cur()
suppressMessages(estimateFwhm(modality = "TIMS-TOF", plot = TRUE, savePlot = svgPath))
stopifnot(grDevices::dev.cur() == device, file.info(svgPath)$size > 0)
estimateFwhm(spectrum, plot = TRUE, savePlot = svgPath)
stopifnot(grDevices::dev.cur() == device, file.info(svgPath)$size > 0)
grDevices::dev.off()
unlink(c(pdfPath, svgPath))
