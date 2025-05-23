% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/useDefaultFwhm.R
\name{useDefaultFwhm}
\alias{useDefaultFwhm}
\title{Use pre-computed default FWHM models}
\usage{
useDefaultFwhm(userMassRes, modality = "FTICR", constant = FALSE)
}
\arguments{
\item{userMassRes:}{a named vector of the form `c(massRes= , mass= )` that specify user-provided mass resolution
at a given m/z. This will used to adjust the height of the FWHM curve depending on the chosen `modality`.}

\item{modality:}{a character specifying the type of the MSI measurement device. Curretly only
`c("FTICR", "TOF-reflector", "TIMS-TOF")` are supported.}

\item{constant:}{a logical, when set to `TRUE`, the fwhm value computed from `userMassRes` is used for the entire
mass range.}
}
\value{
Returns an S3 object 'fwhm' containing a linear interplotor function 'fwhmInterpolator' in addition to
m/z values of the peaks and their corresponding fwhm values.
}
\description{
-- this is still under development --
-- currently only `constant` fwhm values could be provided --
}
\details{
This creates a linear interpolator function that approximats fwhm as a function
of m/z for the given MSI modality.


The defaults are pre-computed and might not represent the true FWHM(m/z) curve for a given dataset.
The user has the option to adjust the hight of the curve if the mass resolution value at a given m/z is provided.
Alternatively, the user could simply provide a single FWHM value which will be used for the entire
mass range. Note that this is not advised but provided (with a warning) for the user's convenience
and for testing.
}
