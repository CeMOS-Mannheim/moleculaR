#' Use default FWHM models
#'
#' Approximate FWHM as a function of m/z when measured peak widths are unavailable.
#'
#' @param userMassRes A numeric vector of the form `c(massRes = , mass = )`,
#'   specifying resolving power and its reference m/z, respectively. Positional
#'   vectors are also accepted. If omitted, the defaults at m/z 249 are 123000
#'   for FTICR, 40000 for TIMS-TOF, and 10000 for TOF-reflector.
#' @param modality Instrument type: `"FTICR"`, `"TOF-reflector"`, or
#'   `"TIMS-TOF"`. The aliases `"TOF"` and `"timsTOF"` are also accepted.
#' @param constant A logical. If `TRUE`, use the width computed from
#'   `userMassRes` across the entire mass range, including for FTICR.
#'
#' @details
#' The FTICR model uses the nonlinear interpolated curve in the bundled
#' `processed-example-Data` dataset. The example instrument's resolving power
#' was 123000 at m/z 249. The curve is multiplied by a positive scale factor
#' so that its width at the reference mass equals `mass / massRes`.
#' This preserves its nonlinear shape without introducing negative widths,
#' which could result from an additive vertical shift. The returned
#' `fwhmVals` contain model predictions rather than the original measured widths.
#'
#' TOF-reflector (e.g. Rapiflex) and TIMS-TOF assume approximately constant
#' resolving power, so FWHM increases linearly with m/z: `FWHM(mz) = mz / massRes`.
#' Their defaults are 10000 and 40000, respectively. The supplied reference
#' mass does not affect this model unless `constant = TRUE`, which explicitly
#' overrides the instrument model with a constant width of `mass / massRes`.
#'
#' These are approximate models. Measured widths from [estimateFwhm()] are
#' preferable when available. The FTICR interpolator retains the bundled
#' curve's endpoint behavior outside its measured mass range.
#'
#' @return An S3 `fwhm` object containing `fwhmInterpolator`, reference
#'   `peaks`, and corresponding model-predicted `fwhmVals` in Da.
#'
#' @examples
#' fticr <- useDefaultFwhm()
#' getFwhm(fticr, 249) # 249 / 123000
#' custom <- useDefaultFwhm(c(massRes = 150000, mass = 400))
#' getFwhm(custom, 400) # 400 / 150000
#' tof <- useDefaultFwhm(modality = "TOF")
#' timstof <- useDefaultFwhm(modality = "timsTOF")
#' getFwhm(timstof, c(249, 800)) # c(249, 800) / 40000
#'
#' @export
useDefaultFwhm <- function(userMassRes = NULL, modality = "FTICR", constant = FALSE) {
      if(!is.character(modality) || length(modality) != 1L || is.na(modality)){
            stop("modality must specify a single instrument type.")
      }
      modality <- switch(modality,
                         "TOF" = "TOF-reflector",
                         "timsTOF" = "TIMS-TOF",
                         modality)
      modality <- match.arg(modality, c("FTICR", "TOF-reflector", "TIMS-TOF"))

      if(!is.logical(constant) || length(constant) != 1L || is.na(constant)){
            stop("constant must be TRUE or FALSE.")
      }
      if(is.null(userMassRes)){
            massRes <- switch(modality, "FTICR" = 123000,
                              "TIMS-TOF" = 40000, "TOF-reflector" = 10000)
            userMassRes <- c(massRes = massRes, mass = 249)
      }
      if(!is.numeric(userMassRes) || length(userMassRes) != 2L ||
         any(!is.finite(userMassRes)) || any(userMassRes <= 0)){
            stop("userMassRes must contain two positive finite numbers: resolving power and reference m/z.")
      }
      if(!is.null(names(userMassRes))){
            if(!setequal(names(userMassRes), c("massRes", "mass"))){
                  stop("Named userMassRes must have names 'massRes' and 'mass'.")
            }
            userMassRes <- userMassRes[c("massRes", "mass")]
      }
      referenceMass <- unname(userMassRes[2])
      width <- unname(referenceMass / userMassRes[1])
      if(!is.finite(width) || width <= 0){
            stop("userMassRes produces a non-positive or non-finite FWHM.")
      }
      if(constant){
            return(fwhm(fwhmInterpolator = approxfun(
                              x = referenceMass, y = width, method = "constant", rule = 2),
                        peaks = referenceMass, fwhmVals = width))
      }
      if(modality != "FTICR"){
            resolvingPower <- unname(userMassRes[1])
            return(fwhm(fwhmInterpolator = function(mzVals) mzVals / resolvingPower,
                        peaks = referenceMass, fwhmVals = width))
      }

      # Load into a private environment to avoid changing the user's workspace.
      exampleData <- new.env(parent = baseenv())
      utils::data(list = "processed-example-Data", package = "moleculaR",
                  envir = exampleData)
      reference <- exampleData$fwhmObj
      rm(exampleData)
      referenceWidth <- getFwhm(reference, referenceMass)
      if(length(referenceWidth) != 1L || !is.finite(referenceWidth) || referenceWidth <= 0){
            stop("The bundled FTICR model has no positive finite width at the reference m/z.")
      }
      scaleFactor <- width / referenceWidth
      if(!is.finite(scaleFactor) || scaleFactor <= 0){
            stop("userMassRes produces an invalid FTICR scale factor.")
      }
      interpolator <- reference$fwhmInterpolator
      scaledInterpolator <- function(mzVals) interpolator(mzVals) * scaleFactor
      peaks <- sort(unique(c(reference$peaks, referenceMass)))
      fwhm(fwhmInterpolator = scaledInterpolator,
           peaks = peaks, fwhmVals = scaledInterpolator(peaks))
}
