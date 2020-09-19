rasterSimulation <-
function(rast, variogramModel) {
\t
\tnscore = function(x)
\t\t{
\t\t\t# Function created by Ashton Shortridge (May/June 2008).\t
   \t\t\t# Takes a vector of values x and calculates their normal scores. Returns a list with the scores
   \t\t\t# and an ordered table of original values and scores, which is useful as a back-transform table.
  \t\t\tnscore = qqnorm(x, plot.it = F)$x  # normal score 
   \t\t\ttrn.table = data.frame(x = sort(x), nscore = sort(nscore))
   \t\t\treturn(list(nscore=nscore, trn.table = trn.table))
\t\t}
\tnscoreBack = function(scores, nscore)
\t\t{
\t\t\t# Function created by Ashton Shortridge (May/June 2008).\t\t
   \t\t\t# Given a vector of normal scores and a normal score object (from nscore), the function returns
   \t\t\t# a vector of back-transformed values.
\t\t\tmin.x = min(nscore$trn.table$x, na.rm=T)
   \t\t\tmax.x = max(nscore$trn.table$x, na.rm=T)
   \t\t\tmin.sc = min(scores, na.rm=T)
   \t\t\tmax.sc = max(scores, na.rm=T)
   \t\t\tx = c(min.x, nscore$trn.table$x, max.x)
   \t\t\tnsc = c(min.sc, nscore$trn.table$nscore, max.sc)
\t\t\tback.xf = approxfun(nsc,x) # Develop the back transform function
   \t\t\tval = back.xf(scores)
\t\t\treturn(val)
\t\t}\t
\tnScoreTransformation = TRUE
\tquantileTransformation = FALSE
\tif (nScoreTransformation == TRUE)
\t\t{
\t\t\trastScore = nscore(rast[])
\t\t\tenvVariableNscore = rast
\t\t\tenvVariableNscore[] = rastScore$nscore
\t\t\tmask = envVariableNscore*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
\t\t\tsimulation = predict(variogramModel, newdata=mask.grid, nsim=1)
\t\t\tsimRaster = raster(simulation)*mask
\t\t\tsimRaster = nscoreBack(simRaster[],rastScore)
\t\t}\t\t\t\t\t\t\t\t
\tif (quantileTransformation == TRUE)
\t\t{
\t\t\tmask = rast*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
\t\t\tsimulation = predict(variogramModel, newdata=mask.grid, nsim=1)
\t\t\tsimRaster = raster(simulation)*mask
\t\t\tenvVariableData = as.data.frame(rast)
\t\t\tnames(envVariableData) = "value"
\t\t\tenvVariableData = envVariableData[with(envVariableData, order(value)),]
\t\t\tsimRasterData = as.data.frame(simRaster)
\t\t\tnames(simRasterData) = "value"
\t\t\tindexes = 1:dim(simRasterData)[1]
\t\t\tsimRasterData$indexes = indexes
\t\t\tsimRasterData = simRasterData[with(simRasterData, order(value)),]
\t\t\tsimRasterData$value = envVariableData
\t\t\tsimRasterData = simRasterData[with(simRasterData, order(indexes)),]
\t\t\tsimRaster = simRasterData[,"value"]
\t\t}
\trast[] = simRaster\t
\treturn(rast)\t
}
