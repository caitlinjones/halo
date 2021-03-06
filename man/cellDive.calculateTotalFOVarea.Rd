% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{cellDive.calculateTotalFOVarea}
\alias{cellDive.calculateTotalFOVarea}
\title{Calculate area of each entire FOV}
\usage{
cellDive.calculateTotalFOVarea(allDat, pp, metaFiles, updatedExclusions,
  allHaloAnnotations = NULL, fovAreaPrefix = "")
}
\arguments{
\item{allDat}{tibble of all Halo data}

\item{pp}{list of all project parameters}

\item{metaFiles}{vector of meta files}

\item{updatedExclusions}{logical indicating whether exclusions were updated during current run of pipeline}

\item{allHaloAnnotations}{list of parsed Halo annotations organized by sample and then FOV}

\item{fovAreaPrefix}{prefix added to 'fov_area_file' and/or 'fov_area_dir' in names of pp list (leave "" if no prefix)}
}
\value{
list of three items:  
          (1) dat = tibble of three columns: Sample, SPOT, Area 
          (2) updated = logical indicating whether a new area table was generated (TRUE) or an old one was loaded (FALSE) 
          (3) pp = updated list of project parameters
}
\description{
Generate a table of three columns: Sample, SPOT, Area
}
