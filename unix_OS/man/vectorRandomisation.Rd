\name{vectorRandomisation}
\alias{vectorRandomisation}
\title{
A function to randomise the end position of a set of movement vectors
}
\description{
This function randomises the position of a set of movement vectors by rotating them around their start location. The randomisation (i.e. rotation) is performed on a template raster and while avoiding that new start and/or new end position fall on a cell with a NODATA ("NA") value. During the randomisation, vectors length (i.e. euclidian distance between start and end locations) is not modified.
}
\usage{
vectorRandomisation(rast, fromCoor, toCoor, showingPlots = FALSE)
}
\arguments{
  \item{rast}{ template raster on which the position of vectors will be randomised. }
  \item{fromCoor}{ matrix containing the coordinates of the start position of each movement vector (one line per vector, 1st column with the longitude and the 2nd column with the latitude). }
  \item{toCoor}{ matrix containing the coordinates of the end position of each movement vector (one line per vector, 1st column with the longitude and the 2nd column with the latitude). }
  \item{showingPlots}{ boolean variable specifying if the plot displaying the randomised vectors has to be displayed or not. }
}
\details{
This function can be used to randomise the movement vectors obtained, e.g., by CMR (capture-mark-recapture) or GPS data.
}
\value{
The function returns a matrix of randomised vector coordinates with one line per vector: new starting longitude, new starting latitude,  new ending longitude and new ending latitude.
}
\references{
Dellicour S, Rose R, Pybus OG (2016a). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016b). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Prunier JG, Piry S, Eloy MC, Bertouille S, Licoppe A, Frantz AC, Flamand MC (2019). Landscape genetic analyses of Cervus elaphus and Sus scrofa: comparative study and analytical developments. \emph{Heredity} 123: 228-241. \cr \cr
Dellicour S, Desmecht D, Paternostre J, Malengreaux C, Licoppe A, Gilbert M, Linden A (2020). Unravelling the dispersal dynamics and drivers of the African swine fever outbreak in Belgium. \emph{Journal of Applied Ecology} 57: 1619-1629.
}
\author{
Simon Dellicour
}
