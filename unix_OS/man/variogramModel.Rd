\name{variogramModel}
\alias{variogramModel}
\title{
A function to estimate variogram model(s)
}
\description{
Using the functions "vgm" and "fit.variogram" of the "gstat" package (Pebesma 2004) to estimate the best variogram model(s) for one or several raster files.
}
\usage{
variogramModel(envVariables)
}
\arguments{
  \item{envVariables}{ list of environmental rasters for which a variogram model has to be estimated. }
}
\value{
The function returns a list of variogram models, one per considered raster file.
}
\author{
Simon Dellicour
}
\note{
This function involves graphics generation and thus cannot be launched from a terminal window. It has to be called from a R interface.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Pebesma EJ (2004). Multivariable geostatistics in S: the gstat package. \emph{Computers & Geosciences} 30: 683-691.
}
