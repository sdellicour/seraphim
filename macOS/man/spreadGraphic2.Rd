\name{spreadGraphic2}
\alias{spreadGraphic2}
\title{
A function to generate a collection of polygons representing the credible intervals of internal node positions.
}
\description{
This function generates a list of spatial polygon dataframes, with one dataframe per considered time slice. Each time slice layer is build by estimating the HPD (highest posterior density) region based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice. See also the software SPREAD of Bielejec \emph{et al}. (2011) for a very similar approach.
}
\usage{
spreadGraphic(localTreesDirectory,
              nberOfExtractionFiles,
              prob = 0.95, 
              startDatum,
              precision = 1)
}
\arguments{
  \item{localTreesDirectory}{ directory where the extracted matrices are saved (see the "treeExctractions" function). }
  \item{nberOfExtractionFiles}{ number of files created with the "treeExtractions" function. }
  \item{prob}{ probability corresponding to the HPD (highest posterior density) regions. }
  \item{startDatum}{ number in a decimal format (e.g. 2007.4) defining the beginning of the dispersal history. }
  \item{precision}{ number, time intervale that will be used to define the successive time slices. }
}
\value{
The function returns a list of spatial polygons dataframe objects, with one object per time slice.
}
\references{
Bielejec F, Rambaut A, Suchard MA, Lemey P (2011). SPREAD: Spatial phylogenetic reconstruction of evolutionary dynamics. \emph{Bioinformatics} 27: 2910-2912. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206.
}
\author{
Simon Dellicour
}
\seealso{
\url{http://www.kuleuven.be/aidslab/phylogeography/SPREAD.html}
}

