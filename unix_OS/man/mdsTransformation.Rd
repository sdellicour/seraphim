\name{mdsTransformation}
\alias{mdsTransformation}
\title{
A function to transform sampling coordinates according to multidimensional scaling (MDS) transformations based on pairwise environmental distances
}
\description{
This function transforms sampling coordinates according to multidimensional scaling (MDS) analyses based on pairwise environmental distances computed among sampling points. A MDS analysis is performed for each specified environmental raster, and original sampling coordinates are replaced by new coordinates found on two main MDS dimensions. Original sampling coordinates have to be contained in sequence names available in a fasta file or in a data frame containing these coordinates (see below for further details). 
}
\usage{
mdsTransformation(input,
                  envVariables,
                  pathModel = 2,
                  resistances = list(),
                  avgResistances = list(),
                  fourCells = FALSE,
                  outputName = "",
                  OS = "Unix")
}
\arguments{
  \item{input}{ list of sequences in a fasta format or data frame containing the sampling coordinates to transform. In the first case, the fasta object has to be obtained by reading a fasta file withe the "read.dna" function of the "ape" package: read.dna ("file_name.fasta", format="fasta"). Important requirement: in the fasta file, the sequence names have to contain the following details in a specific order and separated by underscores: (i) sequence ID (e.g., GenBank accession number), (ii) sampling date in a decimal format (e.g., 2015.47), (iii) latitude coordinate and (iv) longitude coordinate in a decimal format. Here is an example: ">hOH10_1997.2_41.053_-80.706". See the example files for further details. In the second case, the data frame has to contain only four columns: (i) the genomic sequence IDs, (ii) their collection dates, and their sampling (iii) latitude and (iv) longitude coordinates to transform. }
  \item{envVariables}{ list of environmental rasters that will be used to compute the pairwise environmental distances on which MDS analyses will be based. }
  \item{pathModel}{ path model to use for computing environmental distances: "1" (straight-line path model), "2" (least-cost path model, Dijkstra 1959) or "3" (Circuitscape path model, McRae 2006). }
  \item{resistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the rasters of the "envVariables" list have to be treated as resistance or conductance factors (only relevant for the least-cost and Circuitscape path models). }
  \item{avgResistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). This option is only relevant for the Circuitscape path model. }
  \item{fourCells}{ boolean variable specifying if a given raster cell should be connected to either its four first-order (TRUE) or eight first- and second-order (FALSE). This option is only relevant for the least-cost and Circuitscape path models. }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
}
\value{
For each environmental raster that is provided, the function first computes pairwise environmental distances among all sampling coordinates. These pairwise environmental distances are then used to perform MDS analyses and tranform the original geographic coordinates by taking new values on the two main MDS dimensions. In practice, new fasta or text files are generated, one per MDS analysis that is performed. In the case of new fasta files, only the sequence names are modified: longitude and latitude coordinates are updated according to the MDS transformation and thus indicate the new position of each sampled sequence in the transformed space. If text files are generated, these files will then contain the transformed coordinates.
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, \emph{et al.} (\emph{in prep.}). Peaks and valleys of landscape phylogeography. \cr \cr
Dijkstra EW (1959). A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1: 269-271. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561. \cr \cr
Van Etten J (2012). R package gdistance: distances and routes on geographical grids. R package version 1.1.
}
\seealso{
\url{http://cran.r-project.org/web/packages/gdistance/index.html} \cr
\url{http://www.circuitscape.org} \cr
\url{https://github.com/Circuitscape/Circuitscape.jl}
}
