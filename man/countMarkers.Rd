% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_functions.R
\name{countMarkers}
\alias{countMarkers}
\title{Generate a XLSX file of counts for all markers in a given file}
\usage{
countMarkers(markerFile, dataDir, pad = 0, countsXLSXFile = NULL,
  countsRDAFile = NULL, writeXLSXfile = TRUE, saveRDSfile = TRUE,
  runCounts = TRUE, runFracTotal = FALSE, runMedians = TRUE,
  altBases = NULL)
}
\arguments{
\item{markerFile}{File containing list of marker names}

\item{dataDir}{directory of *.rda files containing Halo data}

\item{pad}{amount that will be trimmed from FOV}

\item{countsXLSXFile}{name of XLSX file to write to; if NULL, name will
be automatically generated according to input file names and padding}

\item{countsRDAFile}{name of RDA file to write to; if NULL, RDA file will NOT be written}

\item{writeXLSXfile}{write counts to XLSX file; DEFAULT=TRUE}

\item{saveRDSfile}{save counts table to RDS file; DEFAULT=TRUE}

\item{runCounts}{include counts sheet in output; DEFAULT=TRUE}

\item{runFracTotal}{include fractions of total counts sheet in output; DEFAULT=FALSE}

\item{runMedians}{include medians of counts for each marker in each sample; DEFAULT=TRUE}

\item{altBases}{vector of additional markers for which fractions 
of counts should be calculated; one sheet will be generated
for each element}
}
\description{
Given a marker file and a directory of *.rda files each containing Halo
data for a single sample, generate a XLSX file with multiple sheets of
counts, median counts for each marker in each sample, and fractions of counts
using alternate "bases" (e.g., DAPI)
}
\details{
By default, the file will contain at minimum a sheet of counts and a sheet
containing fractions of total counts. In all sheets, a column represents a marker 
and a row represents a single FOV from a single sample.

If a vector of markers is provided (e.g., c("DAPI","CD3")), another sheet will 
be created for each of them containing the fraction of cells using those markers as 
"bases". IMPORTANT NOTE: Any altBases in addition to DAPI must be counted ONLY IF DAPI
IS POSITIVE AS WELL (e.g., CD3+DAPI+). The names in the output .xlsx file, however
will exclude the DAPI+ to avoid redundancy.
}
