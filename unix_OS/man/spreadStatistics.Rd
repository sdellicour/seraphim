\name{spreadStatistics}
\alias{spreadStatistics}
\title{
A function to compute several dispersion/epidemiological statistics
}
\description{
This function estimates several dispersion/epidemiological statistics based on the spatio-temporal information extracted from phylogenetic trees by the "treeExtractions" function.
}
\usage{
spreadStatistics(localTreesDirectory = ".",
				nberOfExtractionFiles = 1,
				timeSlices = 200, 
				onlyTipBranches = FALSE,
				showingPlots = TRUE, 
				outputName = gsub(" ","_",date()),
				nberOfCores = 1,
				slidingWindow = 1)
}
\arguments{
  \item{localTreesDirectory}{ directory where the extracted matrices are saved (see the "treeExctractions" function). }
  \item{nberOfExtractionFiles}{ number of files created with the "treeExtractions" function. }
  \item{timeSlices}{ number of distinct time slices that will be used to generate the maximal wavefront distance evolution plots. }
  \item{onlyTipBranches}{ boolean variable defining if the statistic estimations have to be performed only while considering the tip branches of the trees. }
  \item{showingPlots}{ boolean variable specifying if the different plots have to be displayed or not. }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{nberOfCores}{ number of available cores to parallelise the computations. }
  \item{slidingWindow}{ sliding window, in units of time, that will be used to generate the dispersal velocity evolution plot (optional). }
}
\value{
The function generates several outputs. The first output is a text file containing the mean branch velocity, the weighted dispersal velocity, the original and weighted diffusion coefficients, the isolation-by-distance (IBD) signal estimated by the Spearman correlation coefficient between the patristic and great-circle geographic distances computed for each pair of tip nodes, as well as branch velocity and diffusion coefficients variation among branches, estimated for each phylogenetic tree. See the related tutorial for the different formulas, as well as Pybus \emph{et al}. (2012) for the formula used to estimate the original diffusion coefficient and Trovao \emph{et al.} (2015) for the formula used to estimate the weighted diffusion coefficient. The function also generates and automatically saves figures displaying the kernel density estimates of mean branch velocity, weighted dispersal velocity and diffusion coefficients parameters (coefficients of variation against mean values). Finally, if the parameter onlyTipBranches = FALSE, the function will also generate and save graphs displaying the evolution of the maximal wavefront distance, as well as the evolution of the dispersal velocity. Two types of wavefront distance are estimated: the maximal \emph{spatial} wavefront distance corresponding to the straight-line distance from to the estimated location of the root, and the maximal \emph{patristic} wavefront distance corresponds to the distance computed as the sum of geographical distances associated with each branch connecting a given node to the root. The maximal \emph{patristic} wavefront distance is here defined as the \emph{patristic} distance from the root to the node associated with the highest \emph{spatial} distance from the root location at a given point in time. On these graphs, grey area corresponds to the 95\% credible region of the estimated wavefront position. Note that geographic distances used to estimate the different statistics are great circle distances computed with the function "rdist.earth" from the R package "fields" (Nychka \emph{et al}. 2014).
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Nychka D, Furrer R, Sain S (2014). fields: Tools for spatial data. R package version 7.1. \cr \cr
Pybus OG, Suchard MA, Lemey P, Bernardin FJ, Rambaut A, Crawford FW, Gray RR, Arinaminpathy N, Stramer SL, Busch MP, Delwart EL (2012). Unifying the spatial epidemiology and molecular evolution of emerging epidemics. \emph{PNAS} 109: 15066-15071. \cr \cr
Trovao NS, Suchard MA, Baele G, Gilbert M, Lemey P (2015). Bayesian inference reveals host-specific contributions to the epidemic expansion of Influenza A H5N1. \emph{Molecular Biology and Evolution} 32 (12): 3264-3275.
}
\author{
Simon Dellicour
}
