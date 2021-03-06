% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_meta.R
\name{getTemplateCellTypeConfig}
\alias{getTemplateCellTypeConfig}
\title{Generate a template marker configuration based on meta data}
\usage{
getTemplateCellTypeConfig(allCellTypes = NULL, cellTypesFile = NULL,
  writeYAML = TRUE, outDir = getwd())
}
\arguments{
\item{allCellTypes}{flattened table containing all cell type file found in meta data files; if NULL,
user must provide cellTypesFile to be parsed; default=NULL}

\item{cellTypesFile}{cell types file in XLSX format structured according to meta data documentation rules;
if NULL, user must provide a pre-generated table consisting of the data found
in meta data file}

\item{writeYAML}{logical indicating whether to write configuration to file; default=TRUE}

\item{outDir}{directory where YAML file should be written; default=getwd()}
}
\value{
all marker information and default marker configuration values
}
\description{
Given either a table of cell types or a cell types file to be parsed, generate in
list format a configuration of marker combinations to be used for downstream 
analysis. The purpose of this template is to provide the beginnings of any and all more complex
configuration needed to calculate statistics and create plots of halo data
}
