\name{mdsTransformation}
\alias{mdsTransformation}
\title{
A function to transform sampling coordinates according to multidimensional scaling (MDS) transformations based on pairwise environmental distances
}
\description{
This function transforms sampling coordinates according to multidimensional scaling (MDS) analyses based on pairwise environmental distances computed among sampling points. A MDS analysis is performed for each specified environmental raster and original sampling coordinates are replaced by new coordinates found on two main MDS dimensions. Original sampling coordinates have to be contained in sequence names available in a fasta file or in tab-delimited text file only containing these coordinates. 
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
  \item{input}{ list of sequences in a fasta format or data frame only containing the sampling coordinates to transform. In the first case, the fasta object has to be obtained by reading a fasta file withe the "read.dna" function of the "ape" package: read.dna ("file_name.fasta", format="fasta"). Important requirement: in the fasta file, the sequence names have to contain the following details in a specific order and separated by underscores: (1) sequence ID (e.g. GenBank accession number), (2) sampling date in a decimal format (e.g. 2015.47), (3) latitude coordinate and (4) longitude coordinate in a decimal format. Here is an example: ">hOH10_1997.2_41.053_-80.706". See the example files for further details. In the second case, the text file has to contain only two columns, one with the latitude and one with the longitude coordinates of each sampled sequence. }
  \item{envVariables}{ list of environmental rasters to that will be used to compute pairwise environmental distances on which MDS analyses will be based. }
  \item{pathModel}{ path taken model to use for computing environmental distances: "1" (straight-line path model), "2" (least-cost path model, Dijkstra 1959) or "3" (Circuitscape path model, McRae 2006). }
  \item{resistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the rasters of the "envVariables" list have to be treated as resistance or conductance factors (only useful for the least-cost and Circuitscape path models). }
  \item{avgResistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). This option is only useful for the Circuitscape path model. }
  \item{fourCells}{ boolean variable specifying if a cell is connected to its four (TRUE) or eight (FALSE) neighbours. This option is only important for the Circuitscape path model. }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
}
\value{
For each environmental raster that is provided, the function first computes pairwise environmental distances among all sampling coordinates. These pairwise environmental distances are then used to perform MDS analyses and tranform the original geographic coordinates by taking new values on the two main MDS dimensions. In practice, new fasta or text files are generated, one per MDS analysis that is performed. In the case of new fasta files, only the sequence names are modified: longitude and latitude coordinates are updated according to the MDS transformation and thus indicate the new position of each sampled sequence in the transformed space. If text files are generated, these files will then only contain the transformed coordinates.
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dijkstra EW (1959). A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1: 269-271. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561. \cr \cr
Van Etten J (2012). R package gdistance: distances and routes on geographical grids. R package version 1.1.
}
\seealso{
\url{http://cran.r-project.org/web/packages/gdistance/index.html} \cr
\url{http://www.circuitscape.org} \cr
}
