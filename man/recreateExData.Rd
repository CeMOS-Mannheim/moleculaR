% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recreateExData.R
\name{recreateExData}
\alias{recreateExData}
\title{Re-create example data}
\usage{
recreateExData(
  pathToImzml,
  pathToSingleSpctr = NULL,
  pathToSldb = NULL,
  pathToMtspc = NULL,
  shrink = TRUE,
  saveTo = getwd()
)
}
\arguments{
\item{pathToImzml:}{path to imzML file.}

\item{pathToSingleSpctr:}{path to the single spectrum tsv file.}

\item{pathToSldb:}{path to the swisslipids database tsv file.}

\item{pathToMtspc:}{path to the Metaspace annotation csv file.}

\item{shrink:}{whether to shrink the msData to only retain the identified m/z in
the Metaspace annotation file.}

\item{saveTo:}{where to save the resulting example data.}
}
\value{
Saves a RData file of the example data.
}
\description{
This function is for internal use only. It re-creates the example data used within the package.
}
\keyword{internal}
