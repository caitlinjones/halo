% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{makeDensityHeatmap}
\alias{makeDensityHeatmap}
\title{Make annotated heatmap of cell type densities}
\usage{
makeDensityHeatmap(den, annot = NULL, annotColors = NULL,
  msToSummarize = NULL, ctConfig = NULL, areas = NULL,
  separateLegend = FALSE, bySample = TRUE, byFOV = TRUE,
  clusterFOVs = TRUE, clusterSamples = TRUE)
}
\arguments{
\item{den}{table containing total FOV cell type densities with basic FOV/Sample annotations}

\item{annot}{table of additional FOV annotation to include in annot tracks}

\item{annotColors}{list of annotation legend colors, where keys match column names in annot, and
each value is a vector of colors named by the respective annotation values}

\item{msToSummarize}{a vector of marker set names to collapse into one "cell type"; if given, must 
provide a parsed cell type configuration and FOV area table; any mutually-exclusive 
(stackable) cell types belonging to a marker set in this vector will be converted to 
a single density value for each FOV.}

\item{ctConfig}{if summarizing any marker sets, must provide a cell type config returned by
getMarkerConfig()}

\item{areas}{if summarizing any marker sets, must provide table of FOV areas in order to calculate
new density.}
}
\description{
Make annotated heatmap of cell type densities
}
