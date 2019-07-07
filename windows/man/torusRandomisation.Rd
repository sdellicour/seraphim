\name{torusRandomisation}
\alias{torusRandomisation}
\title{
A function to perform the torus randomisation of an original raster
}
\description{
This function generates a randomised raster by performing random torus translation and reflection of an original raster (vertical and horizontal reflections both having a probability of occurrence equals to 0.5). Torus translation and reflection randomisation can be used to break down the association between variables while keeping their respective spatial autocorrelation pattern intact, as far as possible (Harms \emph{et al}. 2001, Hardy 2009).
}
\usage{
torusRandomisation(rast)
}
\arguments{
  \item{rast}{ raster to randomise with a torus translation and reflection. }
}
\details{
This function can be called by the "spreadFactor" function.
}
\value{
The function returns a raster object with cell values randomized by a torus translation and reflection.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Harms KE, Condit R, Hubbell SP, Foster RB (2001). Habitat associations of trees and shrubs in a 50-ha neotropical forest plot. \emph{Journal of Ecology} 89: 947-959. \cr \cr
Hardy OJ (2009). TOROCOR: a program to assess the association between spatially autocorrelated variables using a torus-translation test on multiple grids.
}
\author{
Simon Dellicour
}
