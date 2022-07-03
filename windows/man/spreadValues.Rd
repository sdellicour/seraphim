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
				showingPlots = TRUE)
}
\arguments{
  \item{localTreesDirectory}{ directory where the extracted matrices are saved (see the "treeExctractions" function). }
  \item{nberOfExtractionFiles}{ number of files created with the "treeExtractions" function. }
  \item{envVariables}{ list of environmental rasters to analyse. }
  \item{startTime}{ starting time of the plot (x-axis). }
  \item{endTime}{ endingting time of the plot (x-axis). }
  \item{timeSlices}{ number of distinct time slices that will be used to generate the maximal wavefront distance evolution plots. }
  \item{slidingWindow}{ sliding window, in units of time, that will be used to generate the environmental skyline-like plot. }
  \item{showingPlots}{ boolean variable specifying if the different plots have to be displayed or not. }
  \item{outputName}{ name (prefix) to give to the different output files. }
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. }
\author{
Simon Dellicour
}
