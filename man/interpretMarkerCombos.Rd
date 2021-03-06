% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marker_combo_table.R
\name{interpretMarkerCombos}
\alias{interpretMarkerCombos}
\title{Interpret marker combinations}
\usage{
interpretMarkerCombos(markerComboCts, allCellTypes, outFile = NULL)
}
\arguments{
\item{outFile}{if given, markerComboCounts table with interpretations
added will be printed to this file}

\item{markerComboCounts}{table generated by markerComboCounts()}

\item{cellTypesFile}{*.xlsx file defining cell types based on known marker
combination (format in docs)}

\item{unreasonableCombosFile}{*.xlsx file containing a matrix indicating both known 
unreasonable pairwise marker combinations and also
"unknown" pairwise marker combinations}
}
\description{
Given a file of cell types defined by certain marker combinations (format
in docs), assign cell types to each combination in a marker combination counts
file
}
