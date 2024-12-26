\name{isolationByResistance}
\alias{isolationByResistance}
\title{
A function to conduct "isolation-by-resistance" analyses, i.e. to investigate the impact of environmental factors on the isolation-by-distance pattern
}
\description{
This function estimates and tests the correlation between the log-transformed environmental distance(s) and patristic distance associated with each pair of tip nodes in posterior trees inferred through Bayesian continuous phylogeographic inference. While the continuous phylogeographic reconstruction is not used to compute the correlation metric, it is however use to generate a null dispersal model exploited to estimate its statistical support. The function performs three analytical steps: (i) for each environmental factor, the computation of the environmental distances associated with each pair of tip nodes (using the "straight-line path", "least-cost path" or "Circuitscape path" model, see the details below), (ii) the estimation of the correlation between the patristic distances and the corresponding log-transformed environmental distances using a linear regression approach, and (iii) an optional randomisation step to estimate the statistical support of the correlation metric (\emph{Q}). To generate null models, several randomisation procedures of the branches position are implemented (see below for further detail). 
}
\usage{
isolationByResistance(localTreesDirectory = "",
					  nberOfExtractionFiles = 1,
					  envVariables=list(),
					  pathModel = 2,
					  resistances = list(),
					  avgResistances = list(),
					  fourCells = FALSE,
					  nberOfRandomisations = 0,
					  randomProcedure = 3,
					  outputName = "",
					  showingPlots = FALSE,
					  nberOfCores = 1,
					  OS = "Unix",
					  juliaCSImplementation = FALSE)
}
\arguments{
  \item{localTreesDirectory}{ directory where the tree extraction matrices are saved (see the "treeExtractions" or "postTreeExtractions" function). }
  \item{nberOfExtractionFiles}{ number of annotated phylogenetic trees from which spatio-temporal information has been extracted. }
  \item{envVariables}{ list of environmental rasters to test. }
  \item{pathModel}{ path model to use for computing environmental distances: "1" (straight-line path model), "2" (least-cost path model, Dijkstra 1959) or "3" (Circuitscape path model, McRae 2006). }
  \item{resistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the rasters of the "envVariables" list have to be treated as resistance or conductance factors (only relevant for the least-cost and Circuitscape path models). }
  \item{avgResistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). This option is only relevant for the Circuitscape path model. }
  \item{fourCells}{ boolean variable specifying if a given raster cell should be connected to either its four first-order (TRUE) or eight first- and second-order (FALSE). This option is only relevant for the least-cost and Circuitscape path models. }
  \item{nberOfRandomisations}{ number of randomisation steps to perform in order to test the correlation between the patristic and log-transformed environmental distances. }
  \item{randomProcedure}{ randomisation procedure to use: "1" (using external tree randomisations previously performed with the "treesRandomisation" function), "2" (using external RRW simulations previously performed with the "simulatorRRW1" function), "3" (randomisations of nodes position while maintaining the branches length, the tree topology and the location of the most ancestral node), "4" (randomisations of nodes position while maintaining the branches length and the location of the starting node of each branch), "5" (randomisations of nodes position while maintaining the branches length and the location of the ending node of each branch), or "6" (randomisations of nodes position while maintaining only the branches length). }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{showingPlots}{ boolean variable specifying if a plot has to be displayed during the randomisation procedure. }
  \item{nberOfCores}{ number of available cores to parallelise the computations in R (parallelisation only works on Unix systems). }
  \item{nberOfCores_CS}{ number of available cores to parallelise the computations in Circuitscape (we advise to set this parameter value to one and rather use the "numberOfCores" variable to set a parallelisation process). }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
  \item{juliaCSImplementation}{ boolean variable specifying if the Julia implementation of Circuitscape should be called by the function. By default (FALSE), the Python implementation is called. }
}
\value{
The function generates several output files. If the number of randomisations to perform is set to zero, the function will only generate a text file listing the linear regression results for each phylogenetic tree: the \eqn{\beta}{\beta} regression coefficient and coefficient of determination \eqn{R^2_{env}}{R^2_env} obtained from the univariate linear regression between the patristic distances and the log-transformed environmental distances computed for a given environmental factor, as well as the difference \eqn{Q} between \eqn{R^2_{env}}{R^2_env} and \eqn{R^2_{null}}{R^2_null} (the coefficient of determination obtained from the linear regression between the patristic distances and the log-transformed geographical distances associated with each branch). If the number of randomisations to perform is higher than zero, the function will return a second text file reporting the Bayes factor support approximated for each environmental factor and each randomisation step (similar to the approach described in Dellicour \emph{et al}. 2017; see also the related tutorial for detailed information about the Bayes factor estimations). Note that the geographical distance is computed using the selected environmental distance computation method on a "null" raster with uniform cell values equal to "1".
}
\details{
Three different path models are available to compute the environmental distances associated with each pair of tip nodes: (i) the straight-line path model only considering the raster cells below the straight-line segment connecting each pair of tip nodes. In this model, the environmental distance associated with a given pair of tip nodes for a given environmental variable will be the sum of the values of the raster cells below the straight-line segment connecting both tip nodes. (ii) The least-cost path model (Dijkstra 1959). In this second model, environmental distance associated with a given pair of tip nodes for a given environmental variable will be the sum of the values of the raster cells selected by a shortest-path algorithm while considering the environmental variable as a resistance or a conductance factor. (iii) The Circuitscape model based on circuit theory (McRae 2006). This latter model integrates the contribution of multiple possible pathways to compute the environmental distance (McRae 2006). The environmental distances based on the least-cost path model are here computed with the "costDistance" function available in the "gdistance" R package (van Etten 2012) and the environmental distances based on the Circuitscape path model are computed using the Python or Julia package "circuitscape" (McRae 2006). \cr \cr
Regarding the test to estimate the statistical support of the correlation metric \emph{Q}, several randomisation procedures have been implemented to create null models: (i) using already performed tree randomisations (see the "treesRandomisation" function for further detail), (ii) using external RRW simulations corresponding to each tree (see the "simulatorRRW1" function for further detail), (iii) a randomisation of the branches position while maintaining their length, the tree topology and the location of the most ancestral node, (iv) a randomisation of the branches position while maintaining their length and the location of the starting node of each branch, and (v) a randomisation of the branches position while maintaining only their length. The latter procedure can be useful, e.g., when analysing non phylogenetically informed movement data like capture-mark-recapture movement data.
}
\note{
The computation of Circuitscape environmental distances requires the preliminary installation of the program Circuitscape as a Python or Julia package. Instructions on how to install this Python package are available on the Circuitscape website: \url{https://circuitscape.org/}
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Rose R, Faria NR, Vieira LFP, Bourhy H, Gilbert M, Lemey P, Pybus OG (2017). Using viral gene sequences to compare and explain the heterogeneous spatial dynamics of virus epidemics. \emph{Molecular Biology \& Evolution}, \emph{in press}. \cr \cr
Dijkstra EW (1959). A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1: 269-271. \cr \cr
Dellicour S, Troupin C, Jahanbakhsh F, Salama A, Massoudi S, Moghaddam MK, Baele G, Lemey P, Gholami A, Bourhy H (2019). Using phylogeographic approaches to analyse the dispersal history, velocity, and direction of viral lineages – application to rabies virus spread in Iran. \emph{Molecular Ecology} 28: 4335-4350. \cr \cr
Faria NR, Suchard MA, Rambaut A, Streicker DG, Lemey P (2013). Simultaneously reconstructing viral cross-species transmission history and identifying the underlying constraints. \emph{Philosophical Transactions of the Royal Society Series B} 368: 20120196. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561. \cr \cr
Van Etten J (2012). R package gdistance: distances and routes on geographical grids. R package version 1.1.
}
\seealso{
\url{http://cran.r-project.org/web/packages/gdistance/index.html} \cr
\url{http://www.circuitscape.org} \cr
\url{https://github.com/Circuitscape/Circuitscape.jl}
}
