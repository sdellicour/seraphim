\name{spreadFactors}
\alias{spreadFactors}
\title{
A function to investigate the impact of environmental factors on a spread
}
\description{
This function estimates and tests the correlation between environmental distance(s) and dispersal duration associated with each phylogenetic branch. The function performs three steps: (1) for each environmental factor, computation of the environmental distances associated with each phylogenetic branch (using the "straight-line path", "least-cost path" or "Circuitscape" model, see the details below), (2) estimation of the correlation between the dispersal durations and the corresponding environmental distances using a linear regression approach, and (3) optional randomisations to test the level of significance of the different correlation measures. To generate null models, several randomisation procedures of the branches position are implemented. 
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
              OS = "Unix",
              juliaCSImplementation = FALSE)
}
\arguments{
  \item{localTreesDirectory}{ directory where the extracted matrices are saved (see the "treeExctractions" function). }
  \item{nberOfExtractionFiles}{ number of files created with the "treeExtractions" function. }
  \item{envVariables}{ list of environmental rasters to test. }
  \item{pathModel}{ path taken model to use for computing environmental distances: "1" (straight-line path model), "2" (least-cost path model, Dijkstra 1959) or "3" (Circuitscape path model, McRae 2006). }
  \item{resistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the rasters of the "envVariables" list have to be treated as resistance or conductance factors (only useful for the least-cost and Circuitscape path models). }
  \item{avgResistances}{ list of boolean variables corresponding to the "envVariables" list. These boolean variables specify if the raster cells have to be connected by average resistance (TRUE) or average conductance (FALSE). This option is only useful for the Circuitscape path model. }
  \item{fourCells}{ boolean variable specifying if a cell is connected to its four (TRUE) or eight (FALSE) neighbours. This option is only important for the Circuitscape path model. }
  \item{nberOfRandomisations}{ number of randomisation steps to perform in order to test the correlation between environmental distances and dispersal durations. }
  \item{randomProcedure}{ randomisation procedure to use: "1" (using external tree randomisations previously performed with the "treesRandomisation" function), "2" (using external RRW simulations previously performed with the "simulatorRRW1" function), "3" (randomisations of nodes position while maintaining the branches length, the tree topology and the location of the most ancestral node), "4" (randomisations of nodes position while maintaining the branches length and the location of the starting node of each branch), "5" (randomisations of nodes position while maintaining the branches length and the location of the ending node of each branch), or "6" (randomisations of nodes position while maintaining only the branches length). }
  \item{outputName}{ name (prefix) to give to the different output files. }
  \item{showingPlots}{ boolean variable specifying if the different plots have to be displayed or not. }
  \item{OS}{ operating system on which the function is used ("Unix" or "Windows"). }
  \item{juliaCSImplementation}{ boolean variable specifying if the Julia implementation of Circuitscape should be called by the function. By default (FALSE), the Python implementation is called. }
}
\value{
The function creates several output files. If the number of randomisations to perform is set to zero, the function will only create a text file listing the linear regression results for each phylogenetic tree. If the number of randomisations to perform is higher than zero, the function will return a second text file reporting the Bayes factor values approximated for each environmental factor and each randomisation step (Dellicour \emph{et al}. 2017; see also the related tutorial for detailed information about the Bayes factor estimations). Currently, several correlation measures are computed and used: the determination coefficient \eqn{R^2_{env}}{R^2_env} obtained from the univariate linear regression between the dispersal durations and the environmental distances computed for a given environmental factor, its corresponding \emph{p}-value, the difference \eqn{Q} between \eqn{R^2_{env}}{R^2_env} and \eqn{R^2_{null}}{R^2_null} (the determination coefficient obtained from the linear regression between the dispersal durations and the geographical distances associated with each branch). Note that the geographical distance is computed using the selected environmental distance computation method on a "null" raster with uniform cell values equal to 1.  
}
\details{
Three different path taken models are available to compute the environmental distances associated with each phylogenetic branch: (1) the straight-line path model only considering the raster cells below each phylogenetic branch. In this model, the environmental distance associated with a given branch for a given environmental variable will be the sum of the values of the raster cells below the straight-line segment connecting the starting and ending nodes of the branch. (2) The least-cost path model (Dijkstra 1959). In this second model, environmental distance associated with a given branch for a given environmental variable will be the sum of the values of the raster cells selected by a shortest-path algorithm while considering the environmental variable as a resistance or a conductance factor. (3) The Circuitscape model based on circuit theory (McRae 2006). This latter model integrates the contribution of multiple possible pathways to compute the environmental distance (McRae 2006). The least-cost path is here computed with the "costDistance" function available in the "gdistance" R package (van Etten 2012) and the random-walk environmental distance are computed using the Python package "circuitscape" (McRae 2006). \cr \cr
Regarding the significance tests of the correlation measures, several randomisation procedures have been implemented to create null models: (1) using already performed tree randomisations (see the "treesRandomisation" function for more details), (2) using external RRW simulations corresponding to each tree (see the "simulatorRRW1" function for more details), (3) a randomisation of the branches position while maintaining their length, the tree topology and the location of the most ancestral node, (4) a randomisation of the branches position while maintaining their length and the location of the starting node of each branch, and (5) a randomisation of the branches position while maintaining only their length. The latter procedure can be useful, e.g., when analyzing non phylogenetically informed movement data like capture-mark-recapture movement data.
}
\note{
The computation of Circuitscape environmental distances requires the preliminary installation of the software Circuitscape as a Python package. \cr
Instructions on how to install this Python package are available on Circuitscape website: \cr
\url{http://docs.circuitscape.org/circuitscape_4_0_user_guide.html?&id=gsite}
}
\author{
Simon Dellicour
}
\references{
Dellicour S, Rose R, Pybus OG (2016). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. \emph{BMC Bioinformatics} 17: 82. \cr \cr
Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. \emph{Bioinformatics} 32 (20): 3204-3206. \cr \cr
Dellicour S, Rose R, Faria NR, Vieira LFP, Bourhy H, Gilbert M, Lemey P, Pybus OG (2017). Using viral gene sequences to compare and explain the heterogeneous spatial dynamics of virus epidemics. \emph{Molecular Biology \& Evolution}, \emph{in press}. \cr \cr
Dijkstra EW (1959). A note on two problems in connexion with graphs. \emph{Numerische Mathematik} 1: 269-271. \cr \cr
Faria NR, Suchard MA, Rambaut A, Streicker DG, Lemey P (2013). Simultaneously reconstructing viral cross-species transmission history and identifying the underlying constraints. \emph{Philosophical Transactions of the Royal Society Series B} 368: 20120196. \cr \cr
McRae BH (2006). Isolation by resistance. \emph{Evolution} 60: 1551-1561. \cr \cr
Van Etten J (2012). R package gdistance: distances and routes on geographical grids. R package version 1.1.
}
\seealso{
\url{http://cran.r-project.org/web/packages/gdistance/index.html} \cr
\url{http://www.circuitscape.org} \cr
\url{http://ebe.ulb.ac.be/ebe/Software.html}
}
