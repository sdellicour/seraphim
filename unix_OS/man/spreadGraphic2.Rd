\name{spreadGraphic2}
\alias{spreadGraphic2}
\title{
A function to generate a collection of higher posterior density (HPD) regions computed by defined time slice and representing the uncertainty associated with the Bayesian continuous phylogeographic inference
}
\description{
This function generates a list of spatial polygon data frames (SPDFs), with one SPDF per considered time slice. Each time slice layer is build by estimating the HPD (highest posterior density) region based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice. See also, for instance, the programs SPREAD 4 (Nahata \emph{et al}. 2022) and spread.gl (Li \emph{et al}. 2024) implementing similar approaches.
}
\usage{
spreadGraphic2(localTreesDirectory,
               nberOfExtractionFiles,
               prob = 0.95, 
               startDatum,
               precision = 1)
}
\arguments{
  \item{localTreesDirectory}{ directory where the tree extraction matrices are saved (see the "treeExtractions" or "postTreeExtractions" function). }
  \item{nberOfExtractionFiles}{ number of annotated phylogenetic trees from which spatio-temporal information has been extracted. }
  \item{prob}{ probability corresponding to the HPD (highest posterior density) regions. }
  \item{startDatum}{ number in a decimal format (e.g., 2007.4) defining the beginning of the dispersal history. }
  \item{precision}{ number, time interval that will be used to define the successive time slices. }
}
\value{
The function returns a list of spatial polygons dataframe objects, with one object per time slice.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Li Y, Bollen N, Hong SL, Brusselmans M, Gambaro F, Suchard MA, Rambaut A, Lemey P, Dellicour S, Baele G (2024). Spread.gl: visualising pathogen dispersal in a high-performance browser application. \emph{Bioinformatics} 40: btae721. \cr \cr
Nahata K, Bielejec F, Monetta J, Dellicour S, Rambaut A, Suchard MS, Lemey P (2022). SPREAD 4: online visualization of pathogen phylogeographic reconstructions. \emph{Virus Evolution} 8: veac088.
}
\author{
Simon Dellicour
}
