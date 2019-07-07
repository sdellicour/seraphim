\encoding{UTF-8}
\name{readAnnotatedNexus}
\alias{readAnnotatedNexus}
\title{Read annotated tree files in a NEXUS format}
\description{
This function reads annotated trees from files in a NEXUS format. Except for the annotations, this functions mimic ape's function \code{\link[ape]{read.nexus}}. Annotations are ordered by edges, i.e. matching the \code{edge.length} slot of a \code{phylo} object.
}
\usage{
readAnnotatedNexus(file, tree.names=NULL)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character, or a double-quoted string; if ‘file = ""’ (the default) then the tree is input on the keyboard, the entry being terminated with a blank line}
  \item{tree.names}{if there are several trees to be read, a vector of mode character that gives names to the individual trees; if ‘NULL’ (the default), the trees are named ‘"tree1"’, ‘"tree2"’, ...}
}
\details{
  See \code{\link[ape]{read.nexus}} in the ape package for a specification of NEXUS formatted tree files. This function additionally extracts BEAST annotations for all branches/nodes in the trees and returns these annotations as lists of lists in the resulting "phylo" objects
}
\value{
An object of class "phylo" with an additional slot called annotations. This slot is a list indexed by the nodes.
}
\author{
  Marc Suchard
}
\references{
Jombart T, Aanensen DM, Baguelin M, Birrell P, Cauchemez S, Camacho A, Colijn C, Collins C, Cori A, Didelot X, Fraser C, Frost S, Hens N, Hugues J, Höhle M, Opatowski L, Rambaut A, Ratmann O, Soubeyrand S, Suchard MA, Wallinga J, Ypma R, Ferguson N (2014). OutbreakTools: A new platform for disease outbreak analysis using the R software. \emph{Epidemics} 7: 28-34.
}
