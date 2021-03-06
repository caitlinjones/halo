% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{splitByFov}
\alias{splitByFov}
\title{Split data in *.rda file by FOV (SPOT)}
\usage{
splitByFov(dat, outDir)
}
\arguments{
\item{dat}{input tibble/data.frame}

\item{outDir}{output directory}
}
\description{
Take in *.rda file containing data for multiple
FOV and output one *.rda file for each one. Input
tibble must have FOV column named "SPOT". Output
file will be named "FOV_[$SPOT].rda"
}
