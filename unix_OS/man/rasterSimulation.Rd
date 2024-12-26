\name{rasterSimulation}
\alias{rasterSimulation}
\title{
A function to simulate an environmental raster
}
\description{
This function simulates an environmental raster based on a variogram model that can be preliminary estimated with the "variogramModel" function. After the simulation step, a n-score transformation is automatically performed on the simulated raster cell values in order to obtain a distribution of values matching the one observed on the original raster file.
}
\usage{
rasterSimulation(rast, variogramModel)
}
\arguments{
  \item{rast}{ raster template for the raster simulation. }
  \item{variogramModel}{ estimated variogram model for the environmental variable. }
}
\details{
This function was called by a previous version of the "spreadFactor" function.
}
\value{
The function returns a raster simulated under a specific variogram model.
}
\references{
Dellicour S, Rose R, Pybus OG (2016a). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016b). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206.
}
\author{
Simon Dellicour
}
