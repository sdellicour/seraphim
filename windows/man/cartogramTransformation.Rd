\name{cartogramTransformation}
\alias{cartogramTransformation}
\title{
A function to modify sampling coordinates according to cartogram transformation of environmental rasters
}
\description{
This function first performs cartogram transformation of environmental rasters based on their cell values. Cartogram transformations are performed using the method of Gastner & Newman (2004) and while considering each raster cell as a distinct polygon for which an environmental value is assigned (i.e. the original raster cell value). In a second step and for each cartogram transformation that has been performed, the function then modifies sampling geographic coordinates accordingly. These sampling coordinates have to be contained in sequence names available in a fasta file or in tab-delimited text file only containing these coordinates.
}
\usage{
cartogramTransformation(input,
              envVariables,
              resistances = list(),
              outputName = "")           
}
\arguments{
  \item{input}{ list of sequences in a fasta format or data frame only containing the sampling coordinates to transform. In the first case, the fasta object has to be obtained by reading a fasta file withe the "read.dna" function of the "ape" package: read.dna ("file_name.fasta", format="fasta"). Important requirement: in the fasta file, the sequence names have to contain the following details in a specific order and separated by underscores: (1) sequence ID (e.g. GenBank accession number), (2) sampling date in a decimal format (e.g. 2015.47), (3) latitude coordinate and (4) longitude coordinate in a decimal format. Here is an example: ">hOH10_1997.2_41.053_-80.706". See the example files for further details. In the second case, the text file has to contain only two columns, one with the latitude and one with the longitude coordinates of each sampled sequence. }
  \item{envVariables}{ list of environmental rasters to that will be used to perform distinct cartogram transformations. }
  \item{resistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the rasters of the "envVariables" list have to be treated as resistance or conductance factors (only useful for the least-cost and Circuitscape path models). }
  \item{outputName}{ name (prefix) to give to the different output files. }
}
\value{
The function performs one cartogram transformation per provided environmental raster and creates a new fasta/text file per cartogram transformation. In the case of new fasta files, only the sequence names are modified: longitude and latitude coordinates are updated according to the cartogram transformation and thus indicate the new position of each sampled sequence in the transformed space. If text files are generated, these files will then only contain the transformed coordinates.
}
\note{
The use of this function requires the preliminary installation of two additional R packages. These packages can be easily install in R using the "devtools" package: \cr
library(devtools) \cr
system("brew install fftw") \cr
devtools::install_github("omegahat/Rcartogram") \cr
devtools::install_github("chrisbrunsdon/getcartr", subdir="getcartr") \cr
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Gastner MT, Newman MEJ (2004). Diffusion-based method for producing density-equalizing maps. \emph{PNAS} 101: 7499-7504.
}
\seealso{
\url{https://github.com/omegahat/Rcartogram} \cr
\url{https://github.com/chrisbrunsdon/getcartr}
}
