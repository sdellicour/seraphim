\name{treeExtractions}
\alias{treeExtractions}
\title{
A function to extract the spatio-tremporal information embedded in spatiotemporally-referenced phylogenetic trees
}
\description{
This function extracts the spatio-temporal information contained in the phylogenetic trees inferred, e.g., by the continuous phylogeographic method implemented in BEAST (Lemey \emph{et al}. 2010). Depending on the settings of the BEAST analysis, this function can sometimes fail to properly extract information contained in posterior trees. We thus recommend to rather use the latest "postTreeExtractions" function also available in the package.
}
\usage{
treeExtractions(localTreesDirectory,
				allTrees, burnIn,
				randomSampling, 
				nberOfTreesToSample,
				mostRecentSamplingDatum, 
				coordinateAttributeName,
				nberOfCores = 1)
}
\arguments{
  \item{localTreesDirectory}{ name of the (new) directory where the tree extraction matrices will be saved.
  \item{allTrees}{ all the trees in a Nexus format (read with the "scan" function). }
  \item{burnIn}{ number of posterior trees to discard as burn-in. }
  \item{randomSampling}{ boolean variable specifying if the trees have to be randomly sampled in the list. If FALSE, the trees will be sampled at a regular interval (at the largest possible interval). }
  \item{nberOfTreesToSample}{ number of trees to sample. }
  \item{mostRecentSamplingDatum}{ most recent sampling date in a decimal format (e.g., 2007.4). }
  \item{coordinateAttributeName}{ attribute name used to indicate the geographic coordinates within the trees file. }
  \item{nberOfCores}{ number of available cores to parallelise the extractions. }
}
\value{
The function generates "csv" files (one per phylogenetic tree) containing the spatio-tremporal information extracted from annotated phylogenetic trees.
}
\references{
Dellicour S, Rose R, Pybus OG (2016a). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016b). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Lemey P, Rambaut A, Welch JJ, Suchard MA (2010). Phylogeography takes a relaxed random walk in continuous space and time. \emph{Molecular Biology & Evolution} 27: 1877-1885.
}
\author{
Rebecca Rose, Simon Dellicour
}
