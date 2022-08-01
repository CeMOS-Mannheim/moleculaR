#' Use pre-computed default FWHM models
#'
#' -- this is still under development --
#' -- currently only `constant` fwhm values could be provided --
#'
#' This creates a linear interpolator function that approximats fwhm as a function
#' of m/z for the given MSI modality.
#'
#'
#' @param userMassRes: a named vector of the form `c(massRes= , mass= )` that specify user-provided mass resolution
#' at a given m/z. This will used to adjust the height of the FWHM curve depending on the chosen `modality`.
#' @param modality: a character specifying the type of the MSI measurement device. Curretly only
#' `c("FTICR", "TOF-reflector", "TIMS-TOF")` are supported.
#' @param constant:  a logical, when set to `TRUE`, the fwhm value computed from `userMassRes` is used for the entire
#' mass range.
#'
#' @details
#' The defaults are pre-computed and might not represent the true FWHM(m/z) curve for a given dataset.
#' The user has the option to adjust the hight of the curve if the mass resolution value at a given m/z is provided.
#' Alternatively, the user could simply provide a single FWHM value which will be used for the entire
#' mass range. Note that this is not advised but provided (with a warning) for the user's convenience
#' and for testing.
#'
#' @return
#' Returns an S3 object 'fwhm' containing a linear interplotor function 'fwhmInterpolator' in addition to
#' m/z values of the peaks and their corresponding fwhm values.
#'
#' @export
#' @include manualSpatstatImport.R
#'
useDefaultFwhm    <- function(userMassRes, modality = "FTICR", constant = FALSE) {


      # sanity check
      if(length(userMassRes) != 2){
            stop("userMassRes must be a vector of 2 elements; mass resolution and the m/z at which mass resolution was measured.\n")
      }

      # fwhm value at the specified m/z
      fwhmval     <- userMassRes[2] / userMassRes[1]


      if(constant){

            warning("A constant fwhm value is set for the entire mass range, ",
                    "which does not reflect the true fwhm behavious as a function ",
                    "of mass axis. Use this at your own risk. \n")

            return(fwhm(fwhmInterpolator = approxfun(x = userMassRes[2], y = fwhmval, method = "constant", rule = 2),
                        peaks = userMassRes[2],
                        fwhmVals = fwhmval))


      } else {
            stop("This function is still under development. For now only the 'constant=TRUE' is available. \n")
      }


}



