\name{spreadValues}
\alias{spreadValues}
\title{
A function to generate skyline-like plot reporting averaged environmental conditions/values crossed by phylogenetic branches
}
\description{
This function will extract, for each considered time slice, environmental values crossed by the branch segments occurring during this time slice. All the extracted values are then averaged and reported for each sampled tree and per time slice.
}
\usage{
spreadValues(localTreesDirectory = ".",
             nberOfExtractionFiles = 1,
             envVariables = list(),
             startTime,
             endTime,
             timeSlices = 100,
             slidingWindow = 1/12,
             showingPlots = TRUE, 
             outputName = gsub(" ","_",date()),
             nberOfCores = 1)
}
\arguments{
  \item{localTreesDirectory}{ directory where the tree extraction matrices are saved (see the "treeExtractions" or "postTreeExtractions" function). }
  \item{nberOfExtractionFiles}{ number of annotated phylogenetic trees from which spatio-temporal information has been extracted. }
  \item{envVariables}{ list of environmental rasters to analyse. }
  \item{startTime}{ starting time of the plot (x-axis). }
  \item{endTime}{ ending time of the plot (x-axis). }
  \item{timeSlices}{ number of distinct time slices that will be used to generate the maximal wavefront distance evolution plots. }
  \item{slidingWindow}{ sliding window, in units of time, that will be used to generate the environmental skyline-like plot. }
  \item{showingPlots}{ boolean variable specifying if the different plots have to be displayed or not. }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{nberOfCores}{ number of available cores to parallelise the computations (parallelisation will only work on Unix operating systems). }
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Lequime S, Vrancken B, Gill MS, Bastide P, Gangavarapu K, Matteson N, Tan Y, du Plessis L, Fisher AA, Nelson MI, Gilbert M, Suchard MS, Andersen KG, Grubaugh ND, Pybus OG, Lemey P (2020). Epidemiological hypothesis testing using a phylogeographic and phylodynamic framework. \emph{Nature Communications} 11: 5620.
}
\author{
Simon Dellicour
}
