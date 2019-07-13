# List all available and recommended rasters
listAvailableRasters = function() {
    
}

# Download rasters and crop to extent specified. Extent: xmin, xmax, ymin, ymax
downloadRasters = function(ids=c(), extent=c()) {
    
}

# Decrease rasters resolution by a factor 'R'
decreaseResolution = function(envVariables=list(), R=5) {

}

# Scale rasters and return list of 'n' rasters for each scaling k parameter
scaleRasters = function(envVariables=list(), k=c(10,100,1000)) {
    
}

# Download the categorical IGBP land cover and crop to extent specified. Extent: xmin, xmax, ymin, ymax
downloadLandCoverRaster=function(ids=c(), extent=c()) {
    
}

# Generate a distinct continuous raster for each selected land cover variable
prepareLandCoverRasters = function(landCoverRaster, R=5) {
    
}
