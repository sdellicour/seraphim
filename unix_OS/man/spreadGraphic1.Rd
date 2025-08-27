\name{spreadGraphic1}
\alias{spreadGraphic1}
\title{
A function to generate a collection of higher posterior density (HPD) regions retrieved for each internal node from an annotated maximum clade credibility (MCC) tree and representing the uncertainty associated with the continuous phylogeographic inference
}
\description{
This function returns a list of spatial polygon dataframes, with one dataframe per internal node.

See also, for instance, the programs SPREAD 4 (Nahata \emph{et al}. 2022) and spread.gl (Li \emph{et al}. 2024) implementing similar approaches.
}
\usage{
spreadGraphic1(mcc_tre,
               mcc_tab)
}
\arguments{
  \item{mcc_tre}{ MCC tree read/imported by the "readAnnotatedNexus" function of the same package. }
  \item{mcc_tab}{ XXXX. }
}
\details{
XXXX.
}
\value{
XXXX.
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Li Y, Bollen N, Hong SL, Brusselmans M, Gambaro F, Suchard MA, Rambaut A, Lemey P, Dellicour S, Baele G (2024). Spread.gl: visualising pathogen dispersal in a high-performance browser application. \emph{Bioinformatics} 40: btae721. \cr \cr
Nahata K, Bielejec F, Monetta J, Dellicour S, Rambaut A, Suchard MS, Lemey P (2022). SPREAD 4: online visualization of pathogen phylogeographic reconstructions. \emph{Virus Evolution} 8: veac088.
}
\author{
Simon Dellicour
}
