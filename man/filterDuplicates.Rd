% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterDuplicates.R
\name{filterDuplicates}
\alias{filterDuplicates}
\title{Filtersout duplicated analytes}
\usage{
filterDuplicates(spp)
}
\arguments{
\item{spp:}{a spatial point pattern of type 'ppp' with its `marks` being a dataframe with
two columns `idx`and `intensity`.}
}
\description{
Filters out duplicated analytes based on their m/z value. Normally used for the case
of CPPMs i.e. `spp` objects carrying more than one analyte.
}
