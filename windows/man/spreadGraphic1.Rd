\name{spreadGraphic1}
\alias{spreadGraphic1}
\title{
A function to generate a raster file, or a collection of raster files, estimating the credible intervals of internal node positions.
}
\description{
This function generates a raster made by the superimposition of distinct layers corresponding to successive time slices. Each time slice layer is build by estimating the HPD (highest posterior density) region based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice. See also the software SPREAD of Bielejec \emph{et al}. (2011) for a very similar approach. Alternatively, the present function can also generate a collection (list) of raster files each corresponding to a specific time slice.
}
\usage{
spreadGraphic(localTreesDirectory,
              nberOfExtractionFiles,
              rast,
              prob = 0.95, 
              startDatum,
              precision = 1,
              timeLayers = FALSE,
              origin = FALSE)
}
\arguments{
  \item{localTreesDirectory}{ directory where the extracted matrices are saved (see the "treeExctractions" function). }
  \item{nberOfExtractionFiles}{ number of files created with the "treeExtractions" function. }
  \item{rast}{ template raster that will be used to map the dispersal history (raster cell values will not be taken into account). }
  \item{prob}{ probability corresponding to the HPD (highest posterior density) regions. }
  \item{startDatum}{ number in a decimal format (e.g. 2007.4) defining the beginning of the dispersal history. }
  \item{precision}{ number, time intervale that will be used to define the successive time slices. }
  \item{timeLayers}{ boolean variable specifying if the function has to generate a single raster ("timeLayers=FALSE", default) or a list of raster files with one file per considered time slice. }
  \item{origin}{ boolean variable specifying if the function has to rather estimate the HPD region of the spread origin (default value is "FALSE"; see the details). }
}
\details{
If the "origin" boolean variable is set to "TRUE", the function will only estimate the HPD region of the spread origin. In that particular case, the function will only focus on the position of the most ancestral node of each extraction matrix, no mater their associated time. This latter aspect differs from the main option ("origin = FALSE", i.e. estimation of the dispersal history).
}
\value{
The function returns a raster object or a list of raster objects.
}
\references{
Bielejec F, Rambaut A, Suchard MA, Lemey P (2011). SPREAD: Spatial phylogenetic reconstruction of evolutionary dynamics. \emph{Bioinformatics} 27: 2910-2912. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206.
}
\author{
Simon Dellicour, Nuno Faria
}
\seealso{
\url{http://www.kuleuven.be/aidslab/phylogeography/SPREAD.html}
}

