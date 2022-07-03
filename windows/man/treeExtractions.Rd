\name{treeExtractions}
\alias{treeExtractions}
\title{
A function to extract the spatio-tremporal information contained in spatiotemporally-referenced phylogenetic trees
}
\description{
This function extracts the spatio-temporal information contained in the phylogenetic trees inferred e.g. by the continuous phylogeographic method implemented in BEAST (Lemey \emph{et al}. 2010).
}
\usage{
treeExtractions(localTreesDirectory,
                allTrees, burnIn,
                randomSampling, 
                nberOfTreesToSample,
                mostRecentSamplingDatum, 
                coordinateAttributeName)
}
\arguments{
  \item{localTreesDirectory}{ name of the (new) directory where the extracted matrices have to be saved. }
  \item{allTrees}{ all the trees in a Nexus format. }
  \item{burnIn}{ number of trees to discard as burn-in. }
  \item{randomSampling}{ boolean variable specifying if the trees have to be randomly sampled in the list. If FALSE, the trees will be sampled at a regular interval (at the largest possible intervale). }
  \item{nberOfTreesToSample}{ number of trees to sample. }
  \item{mostRecentSamplingDatum}{ most recent datum of sampling in decimal format (e.g. 2007.4). }
  \item{coordinateAttributeName}{ attribute name used to indicate the geographic coordinates within the trees file. }
}
\value{
The function generates "csv" files (one per phylogenetic tree) containing the spatio-tremporal information extracted from phylogenetic trees.
}
\references{
Dellicour S, Rose R, Pybus OG (2016). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Lemey P, Rambaut A, Welch JJ, Suchard MA (2010). Phylogeography takes a relaxed random walk in continuous space and time. \emph{Molecular Biology & Evolution} 27: 1877-1885.
}
\author{
Rebecca Rose, Simon Dellicour
}
