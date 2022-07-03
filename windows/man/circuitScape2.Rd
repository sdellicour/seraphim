\name{circuitScape2}
\alias{circuitScape2}
\title{
A function to launch Circuitscape analyses
}
\description{
The function calls the external Julia library "circuitscape" to compute resistance or conductance environmental distance associated with each phylogenetic branch.
}
\usage{
circuitScape1(envVariable,
              envVariableName,
              resistance = TRUE,
              avgResistance = TRUE,
              fourCells = FALSE,
              fromCoor,
              toCoor,
              OS = "Unix",
              prefix = "",
              ID = "")
}
\arguments{
\item{envVariable}{ raster file on which the resistance or conductance has to be computed. }
  \item{envVariableName}{ name (string) of the raster file. This name can be different from the name saved in the raster object and will be used to name the temporary file generated for the Circuitscape analysis. }
  \item{resistance}{ boolean variable specifying if the raster will be treated as a resistance or conductance factor. }
  \item{avgResistance}{ boolean variable specifying if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). }
  \item{fourCells}{ boolean variable specifying if a cell is connected to its four (TRUE) or eight (FALSE) neighbours. }
  \item{fromCoor}{ matrix containing the coordinates of the starting nodes of each phylogenetic branch coming from one tree (one line per node, 1st column with the longitude and the 2nd column with the latitude). }
  \item{toCoor}{ matrix containing the coordinates of the ending nodes of each phylogenetic branch coming from one tree (one line per node, 1st column with the longitude and the 2nd column with the latitude). }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
  \item{prefix}{ a character string that will be used as a prefix to name the temporary folders generated for Circuitscape. }
  \item{ID}{ a number to identify the Circuitscape operation. }
}
\value{
The function returns a vector containing the resistance or conductance estimated for each phylogenetic branch.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561.
}
\author{
Simon Dellicour
}
\note{
This function requires the preliminary installation of the software Circuitscape as a Pyhton package. Instructions on how to install this Python package are available on Circuitscape website: \cr
\url{http://docs.circuitscape.org/circuitscape_4_0_user_guide.html}. \cr
Also note that fromCoor and toCoor matrices have to display the same dimensions and lines in the same order: coordinates in line i in fromCoor have to be the coordinates of the starting node connected to the ending node for which geographic coordinates are on line i in toCoor. \cr
This function can be used the Python library "circuitscape" to compute resistance or conductance environmental distance associated with any movement vectors for which start and end locations are respectively specified in the "fromCoor" and "toCoor" arguments.
}
