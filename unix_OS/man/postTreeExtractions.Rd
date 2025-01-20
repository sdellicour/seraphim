\name{postTreeExtractions}
\alias{postTreeExtractions}
\title{
A function to extract the spatio-tremporal information embedded in annotated phylogenetic trees sampled from the posterior distribution of a Bayesian continuous phylogeographic inference
}
\description{
This function extracts the spatio-temporal information contained in the annotated phylogenetic trees inferred, e.g., by the continuous phylogeographic method implemented in BEAST (Lemey \emph{et al}. 2010). To use this function the "latitute-longitude" continuous traits partition should have been named "location" while preparing the input XML file for the continuous phylogeographic analysis conducted with BEAST.
}
\usage{
postTreeExtractions(post_tre, mostRecentSamplingDatum)
}
\arguments{
  \item{post_tre}{ annotated phylogenetic tree read/imported by the "readAnnotatedNexus" function of the same package. }
  \item{mostRecentSamplingDatum}{ most recent sampling date in a decimal format (e.g., 2007.4). }
}
\value{
The function returns an extraction table containing the spatio-tremporal information extracted from an annotated phylogenetic tree.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Lemey P, Rambaut A, Welch JJ, Suchard MA (2010). Phylogeography takes a relaxed random walk in continuous space and time. \emph{Molecular Biology & Evolution} 27: 1877-1885.
}
\author{
Simon Dellicour
}
