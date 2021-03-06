% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/melanoma_spatial.R
\name{calculateBandAreas}
\alias{calculateBandAreas}
\title{Calculate area of interface bands, each defined as the collection of points that are located
between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]}
\usage{
calculateBandAreas(allBoundaries, maxG = 5, bb, interfaceBins)
}
\arguments{
\item{allBoundaries}{halo boundaries after removal of contained boundaries (list returned by 
cleanBoundaries())}

\item{maxG}{maximum factor of 10 to use for generating random points on grid [ REWORD??? ]}

\item{bb}{a list containing X0,X1,Y0,Y1 representing the boundary
of the trimmed FOV}

\item{interfaceBins}{vector of distances (integers) that define each band; default = (-20:20)*10}
}
\value{
a dataframe with a row for each interface bin and a column for each area calculation (one for 
        10^maxG and one for 10^maxG-1)
}
\description{
Calculate area of each band, defined as the collection of cells that are located
between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]
**Do this for 10^maxG and for 10^maxG-1
}
