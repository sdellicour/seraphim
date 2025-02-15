\name{treesRandomisation}
\alias{treesRandomisation}
\title{
A function to randomise the position of pmhylogenetic tree branches
}
\description{
This function randomises the position of phylogenetic branches on the study area and is based "extraction" tables previously generated by the "treeExtractions" or "postTreeExtractions" function, which both summarise each phylogenetic branch as a movement vector. Four distinct randomisations have been implemented (see below).
}
\usage{
treesRandomisation(localTreesDirectory = "",
                   nberOfExtractionFiles = 1,
                   envVariables=list(),
                   randomProcedure = 3,
                   nberOfCores = 1)
}
\arguments{
  \item{localTreesDirectory}{ directory where the tree extraction matrices are saved (see the "treeExtractions" or "postTreeExtractions" function). }
  \item{nberOfExtractionFiles}{ number of annotated phylogenetic trees from which spatio-temporal information has been extracted. }
  \item{envVariables}{ list of environmental rasters to consider. These rasters are only use to prevent that randomised node positions fall in non-accessible areas, i.e. raster cells associated with "NA" value (e.g., sea areas). }
  \item{randomProcedure}{ randomisation procedure to use: "3" (randomisations of nodes position while maintaining the branches length, the tree topology and the location of the most ancestral node), "4" (randomisations of nodes position while maintaining the branches length and the location of the starting node of each branch), "5" (randomisations of nodes position while maintaining the branches length and the location of the ending node of each branch), or "6" (randomisations of nodes position while maintaining only the branches length). }
  \item{nberOfCores}{ number of available cores to parallelise the computations in R (parallelisation will only work on Unix operating systems). }
}
\value{
In practice, the function generates a table containing the spatio-tremporal information corresponding to each randomised tree, and which are similar to the tables generated by the "treeExtractions" or "postTreeExtractions" function.  
}
\details{
Several randomisation procedures of tree branches have been implemented: (i) a randomisation of the branches position while maintaining their length, the tree topology and the location of the most ancestral node (procedure n°3), (ii) a randomisation of the branches position while maintaining their length and the location of the starting node of each branch (procedure n°4), and (iii) a randomisation of the branches position while maintaining only their length  (procedure n°5). The latter procedure can be useful, e.g., when analysing non phylogenetically informed movement data like capture-mark-recapture movement data. Note that in each case, branch positions are independently randomised for each tree.
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016a). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Rose R, Pybus OG (2016b). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82.
}
