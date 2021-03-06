% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{cellDive.loadAllData}
\alias{cellDive.loadAllData}
\title{Load data from each rda file into one table}
\usage{
cellDive.loadAllData(projectParams, whichFileSet)
}
\arguments{
\item{projectParams}{key-value pairs of project parameters in list format}

\item{whichFileSet}{name of key (minus "_dir" and "_files" suffixes) of files to 
be loaded (e.g., 'raw_data' or 'data' or 'rethresholded_data')}
}
\value{
tibble containing all data required for downstream analysis
}
\description{
Read each data file, select the minimum columns required for downstream
analyses and store as one table in memory. Columns include: 
    Sample
    SubSample
    SPOT
    UUID
    Marker
    XMin/Max
    YMin/Max
    Value
Data is filtered for cells whose EXCLUDE value == "" and ValueType == "Positive"
}
