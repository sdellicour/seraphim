\name{plotRaster}
\alias{plotRaster}
\title{
A function to plot a raster with a simplified layout
}
\description{
This function generates a raster plot using a simplified layout.
}
\usage{
plotRaster(rast,
           cols = NULL,
           colNA = "gray90",
           add = FALSE,
           new = TRUE, 
           addBox = TRUE,
           addAxes = FALSE,
           addLegend = FALSE,
           legendOnly = FALSE)
}
\arguments{
  \item{rast}{ raster object to plot. }
  \item{cols}{ colour scale to use (by default, it will be the "RdBu" colour scale from the "RColorBrewer" package. }
  \item{colNA}{ colour to use for the "NA" raster cell values (default is "gray90"). }
  \item{add}{ boolean variable specifying if the plot has to be added to the current plot. } 
  \item{new}{ boolean variable specifying if the plot has to be displayed in a new window (only used if "add" is set to "FALSE"). }
  \item{addBox}{ boolean variable specifying if a box has to be added around the raster plot (only used if "add" is set to "FALSE"). }
  \item{addAxes}{ boolean variable specifying if axes displaying minimum and maximum geographic coordinates have to added to the plot (only used if "add" is set to "FALSE"). }
  \item{addLegend}{ boolean variable specifying if a legend has to be added to the plot. }
  \item{legendOnly}{ boolean variable specifying if only the legend has to be added to a plot. }
}
\author{
Simon Dellicour
}
