\name{mccTreeExtractions}
\alias{mccTreeExtractions}
\title{
A function to extract the spatio-tremporal information embedded in the maximum clade credibility (MCC) tree obtained from a continuous phylogeographic inference conducted with the software package BEAST
}
\description{
This function extracts the spatio-temporal information contained in the maximum clade credibility (MCC) phylogenetic trees inferred by the continuous phylogeographic method implemented in BEAST (Lemey \emph{et al}. 2010). To use this function the "latitute-longitude" continuous traits partition should have been named "location" while preparing the input XML file for the continuous phylogeographic analysis conducted with BEAST.
}
}
\usage{
mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)
}
\arguments{
  \item{mcc_tre}{ MCC tree read/imported by the "readAnnotatedNexus" function of the same package. }
  \item{mostRecentSamplingDatum}{ most recent sampling date in a decimal format (e.g., 2007.4). }
}
\value{
The function returns an extraction table containing the spatio-tremporal information extracted from the MCC tree.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Lemey P, Rambaut A, Welch JJ, Suchard MA (2010). Phylogeography takes a relaxed random walk in continuous space and time. \emph{Molecular Biology & Evolution} 27: 1877-1885.
}
\author{
Simon Dellicour
}
