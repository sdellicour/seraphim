\name{spreadGraphic1}
\alias{spreadGraphic1}
\title{
A function to generate a raster file, or a collection of raster files, estimating the credible intervals of internal node positions inferred by a continuous phylogeographic analysis
}
\description{
This function generates a raster made by the superimposition of distinct layers corresponding to successive time slices. Each time slice layer is build by estimating the HPD (highest posterior density) region based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice. See also, for instance, the programs SPREAD 4 (Nahata \emph{et al}. 2022) and spread.gl (Li \emph{et al}. 2024) for a very similar approach. Alternatively, the present function can also generate a collection (list) of raster files each corresponding to a specific time slice.
}
\usage{
spreadGraphic1(localTreesDirectory,
               nberOfExtractionFiles,
               rast,
               prob = 0.95, 
               startDatum,
               precision = 1,
               timeLayers = FALSE,
               nberOfCores = 1, 
               origin = FALSE)
}
\arguments{
  \item{localTreesDirectory}{ directory where the tree extraction matrices are saved (see the "treeExtractions" or "postTreeExtractions" function). }
  \item{nberOfExtractionFiles}{ number of annotated phylogenetic trees from which spatio-temporal information has been extracted. }
  \item{rast}{ template raster that will be used to map the HPD regions (raster cell values will not be taken into account). }
  \item{prob}{ probability corresponding to the HPD regions. }
  \item{startDatum}{ number in a decimal format (e.g., 2007.4) defining the beginning of the dispersal history. }
  \item{precision}{ number, time interval that will be used to define the successive time slices. }
  \item{timeLayers}{ boolean variable specifying if the function has to generate a single raster ("timeLayers=FALSE", default) or a list of raster files with one file per considered time slice. }
  \item{nberOfCores}{ number of available cores to parallelise the estimation of the different time slice layers (parallelisation only works on Unix systems). }
  \item{origin}{ boolean variable specifying if the function has to rather estimate the HPD region of the spread origin (default value is "FALSE"). }
}
\details{
If the "origin" boolean variable is set to "TRUE", the function will only estimate the HPD region of the spread origin. In that particular case, the function will only focus on the position of the most ancestral node of each extraction matrix, no mater their associated time. This latter aspect differs from the main option ("origin = FALSE", i.e. estimation of the dispersal history).
}
\value{
The function returns a raster object or a list of raster objects.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Li Y, Bollen N, Hong SL, Brusselmans M, Gambaro F, Suchard MA, Rambaut A, Lemey P, Dellicour S, Baele G (2024). Spread.gl: visualising pathogen dispersal in a high-performance browser application. \emph{Bioinformatics} 40: btae721. \cr \cr
Nahata K, Bielejec F, Monetta J, Dellicour S, Rambaut A, Suchard MS, Lemey P (2022). SPREAD 4: online visualization of pathogen phylogeographic reconstructions. \emph{Virus Evolution} 8: veac088.
}
\author{
Simon Dellicour, Nuno Faria
}
