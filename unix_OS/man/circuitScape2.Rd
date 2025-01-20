\name{circuitScape2}
\alias{circuitScape2}
\title{
A function to launch Circuitscape analyses with the related Julia package
}
\description{
The function calls the external Julia library "circuitscape" to compute the environmental distance associated with each phylogenetic branch of a spatially annotated tree, e.g., inferred by a continuous phylogeographic analysis.
}
\usage{
circuitScape2(envVariable,
              envVariableName,
              resistance = TRUE,
              avgResistance = TRUE,
              fourCells = FALSE,
              fromCoor,
              toCoor,
              OS = "Unix",
              prefix = "",
              ID = "",
              nberOfCores_CS=1)
}
\arguments{
\item{envVariable}{ raster file on which environmental distances have to be computed. }
  \item{envVariableName}{ name (string) of the raster file. This name can be different from the name saved in the raster object and will be used to name the temporary file generated for the Circuitscape analysis. }
  \item{resistance}{ boolean variable specifying if the raster will be treated as a resistance or conductance factor. }
  \item{avgResistance}{ boolean variable specifying if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). }
  \item{fourCells}{ boolean variable specifying if a given raster cell should be connected to either its four first-order (TRUE) or eight first- and second-order (FALSE). This option is only relevant for the least-cost and Circuitscape path models. }
  \item{fromCoor}{ matrix containing the coordinates of the starting nodes of each phylogenetic branch coming from one tree (one line per node, 1st column with the longitude and the 2nd column with the latitude). }
  \item{toCoor}{ matrix containing the coordinates of the ending nodes of each phylogenetic branch coming from one tree (one line per node, 1st column with the longitude and the 2nd column with the latitude). }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
  \item{prefix}{ a character string that will be used as a prefix to name the temporary folders generated for Circuitscape. }
  \item{ID}{ a number to identify the Circuitscape operation. }
  \item{nberOfCores_CS}{ number of available cores to parallelise the computations in Circuitscape (we advise to set this parameter value to one and rather use the "numberOfCores" variable to set a parallelisation process). }
}
\value{
The function returns a vector containing the environmental distance computed for each phylogenetic branch.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016a). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Rose R, Pybus OG (2016b). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561.
}
\author{
Simon Dellicour
}
\note{
This function requires the preliminary installation of the program Circuitscape as a Julia package. Instructions on how to install this Julia package are available here: \cr
https://docs.circuitscape.org/Circuitscape.jl/latest/}. \cr
Also note that "fromCoor" and "toCoor" matrices have to display the same dimensions and lines in the same order: coordinates in line \emph{i} in "fromCoor" have to be the coordinates of the starting node connected to the ending node for which geographic coordinates are on line \emph{i} in "toCoor"". \cr
This function can be used to compute (resistance or conductance) environmental distance associated with any movement vector for which start and end locations are respectively specified in the "fromCoor" and "toCoor" arguments.
}
