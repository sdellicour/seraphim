\name{spreadFactors}
\alias{spreadFactors}
\title{
A function to investigate the impact of environmental factors on the dispersal dynamic of (viral) lineages
}
\description{
This function estimates and tests the correlation between environmental distance(s) and dispersal duration associated with each phylogenetic branch. The function performs three analytical steps: (i) for each environmental factor, the computation of the environmental distances associated with each phylogenetic branch (using the "straight-line path", "least-cost path" or "Circuitscape path" model, see the details below), (ii) the estimation of the correlation between the dispersal durations and the corresponding environmental distances using a linear regression approach, and (iii) an optional randomisation step to estimate the statistical support of the correlation metrics (\emph{Q}). To generate null models, several randomisation procedures of the branches position are implemented (see below for further detail). 
}
\usage{
spreadFactors(localTreesDirectory = "",
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
  \item{nberOfRandomisations}{ number of randomisation steps to perform in order to test the correlation between environmental distances and dispersal durations. }
  \item{randomProcedure}{ randomisation procedure to use: "1" (using external tree randomisations previously performed with the "treesRandomisation" function), "2" (using external RRW simulations previously performed with the "simulatorRRW1" function), "3" (randomisations of nodes position while maintaining the branch lengths, the tree topology, and the location of the most ancestral node), "4" (randomisations of nodes position while maintaining the branch lengths and the location of the starting node of each branch), "5" (randomisations of node positions while maintaining the branch lengths and the location of the ending node of each branch), or "6" (randomisations of node positions while maintaining only the branch lengths). }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{showingPlots}{ boolean variable specifying if a plot has to be displayed during the randomisation procedure. }
  \item{nberOfCores}{ number of available cores to parallelise the computations in R (parallelisation will only work on Unix operating systems). }
  \item{nberOfCores_CS}{ number of available cores to parallelise the computations in Circuitscape (we advise to set this parameter value to one and rather use the "numberOfCores" variable to set a parallelisation process). }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
  \item{juliaCSImplementation}{ boolean variable specifying if the Julia implementation of Circuitscape should be called by the function. By default (FALSE), the Python implementation is called. }
}
\value{
The function generates several output files. If the number of randomisations to perform is set to zero, the function will only generate a text file listing the results of the two linear regression (LR) analyses conducted for each phylogenetic tree: (i) the univariate linear regression between the dispersal durations (\eqn{t}) and the environmental distances (\eqn{d}) computed for a given environmental factor (for which the output files are labelled "LR1"), and (ii) the univariate linear regression between four times the dispersal durations (\eqn{4t}) and the squared environmental distances (\eqn{d^2}) computed for a given environmental factor (for which the output files are labelled "LR2"). While the first linear regression (LR1) aims to assess the impact of a given environmental factor on the heterogeneity of the lineage dispersal velocity (\eqn{t \sim d}), the second linear regression (LR2) aims to assess the impact of a given environmental factor on the heterogeneity of the diffusion velocity (\eqn{4t \sim d^2}), which is related to the diffusion coefficient identified as a metric more robust to the sampling intensity (i.e. the sampling size) than the lineage dispersal velocity (Dellicour \emph{et al}. 2024). For both linear regressions, the function will thus generate a first text file gathering the following elements: the \eqn{\beta}{\beta} regression coefficients and coefficient of determination \eqn{R^2_{env}}{R^2_env} obtained from each linear regression, as well as the difference \eqn{Q} between \eqn{R^2_{env}}{R^2_env} and \eqn{R^2_{null}}{R^2_null} (the coefficient of determination obtained from the linear regression based on environmental distances computed on a "null" raster with uniform cell values equal to "1", the later distances corresponding to a proxy of the geographic distances). If the number of randomisations to perform is higher than zero, the function will return a additional text files reporting the Bayes factor support approximated for each environmental factor and each randomisation step (Dellicour \emph{et al}. 2017; see also the related tutorial for detailed information about the Bayes factor estimations). Note that if the "pathModel" parameter is set to "0", the "spreadFactors" function will rather conduct association analyses between the lineage dispersal locations and local environmental distances (Dellicour \emph{et al}. 2019) — those analyses and their outcomes are not detailed here but in a dedicated tutorial ("impact on lineage location").
}
\details{
Three different path models are available to compute the environmental distances associated with each phylogenetic branch: (i) the straight-line path model only considering the raster cells below each phylogenetic branch. In this model, the environmental distance associated with a given branch for a given environmental variable will be the sum of the values of the raster cells below the straight-line segment connecting the starting and ending nodes of the branch. (ii) The least-cost path model (Dijkstra 1959). In this second model, environmental distance associated with a given branch for a given environmental variable will be the sum of the values of the raster cells selected by a shortest-path algorithm while considering the environmental variable as a resistance or a conductance factor. (iii) The Circuitscape model based on circuit theory (McRae 2006). This latter model integrates the contribution of multiple possible pathways to compute the environmental distance (McRae 2006). The environmental distances based on the least-cost path model are here computed with the "costDistance" function available in the "gdistance" R package (van Etten 2012) and the environmental distances based on the Circuitscape path model are computed using the Python or Julia package "circuitscape" (McRae 2006). \cr \cr
Regarding the test to estimate the statistical support of the correlation metrics \emph{Q}, several randomisation procedures have been implemented to create null models: (i) using already performed tree randomisations corresponding to each tree (see the "treesRandomisation" function for further detail), (ii) using already performed RRW simulations corresponding to each tree (see the "simulatorRRW1" function for further detail), (iii) a randomisation of the branches position while maintaining their length, the tree topology and the location of the most ancestral node, (iv) a randomisation of the branches position while maintaining their length and the location of the starting node of each branch, and (v) a randomisation of the branches position while maintaining only their length. The latter procedure can be useful, e.g., when analysing non phylogenetically informed movement data like capture-mark-recapture movement data.
}
\note{
The computation of Circuitscape environmental distances requires the preliminary installation of the program Circuitscape as a Python or Julia package. Instructions on how to install this Python package are available on the Circuitscape website:
\url{https://circuitscape.org/}
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016a). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Rose R, Pybus OG (2016b). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria NR, Vieira LFP, Bourhy H, Gilbert M, Lemey P, Pybus OG (2017). Using viral gene sequences to compare and explain the heterogeneous spatial dynamics of virus epidemics. \emph{Molecular Biology \& Evolution} 34: 2563-2571. \cr \cr
Dellicour S, Troupin C, Jahanbakhsh F, Salama A, Massoudi S, Moghaddam MK, Baele G, Lemey P, Gholami A, Bourhy H (2019). Using phylogeographic approaches to analyse the dispersal history, velocity, and direction of viral lineages – application to rabies virus spread in Iran. \emph{Molecular Ecology} 28: 4335-4350. \cr \cr
Dellicour S, Bastide P, Rocu P, Fargette D, Hardy OJ, Suchard MA, Guindon S, Lemey P (2024). How fast are viruses spreading in the wild? \emph{PLoS Biology} 22: e3002914. \cr \cr
Dellicour S, Gámbaro F, Jacquot M, Lequime S, Baele G, Gilbert M, Pybus OG, Suchard MA, Lemey P (2025). \emph{Comparative performance of viral landscape phylogeography approaches. Proceedings of the National Academy of Sciences of the USA} 122: e2506743122. \cr \cr
Dijkstra EW (1959). A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1: 269-271. \cr \cr
Faria NR, Suchard MA, Rambaut A, Streicker DG, Lemey P (2013). Simultaneously reconstructing viral cross-species transmission history and identifying the underlying constraints. \emph{Philosophical Transactions of the Royal Society Series B} 368: 20120196. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561. \cr \cr
Van Etten J (2012). R package gdistance: distances and routes on geographical grids. R package version 1.1.
}
\seealso{
\url{https://beast.community/} \cr
\url{http://cran.r-project.org/web/packages/gdistance/index.html} \cr
\url{http://www.circuitscape.org} \cr
\url{https://github.com/Circuitscape/Circuitscape.jl}
}
