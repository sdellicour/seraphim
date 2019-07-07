rasterSimulation <-
function (rast, variogramModel) 
{
    nscore = function(x) {
        nscore = qqnorm(x, plot.it = F)$x
        trn.table = data.frame(x = sort(x), nscore = sort(nscore))
        return(list(nscore = nscore, trn.table = trn.table))
    }
    nscoreBack = function(scores, nscore) {
        min.x = min(nscore$trn.table$x, na.rm = T)
        max.x = max(nscore$trn.table$x, na.rm = T)
        min.sc = min(scores, na.rm = T)
        max.sc = max(scores, na.rm = T)
        x = c(min.x, nscore$trn.table$x, max.x)
        nsc = c(min.sc, nscore$trn.table$nscore, max.sc)
        back.xf = approxfun(nsc, x)
        val = back.xf(scores)
        return(val)
    }
    nScoreTransformation = TRUE
    quantileTransformation = FALSE
    if (nScoreTransformation == TRUE) {
        rastScore = nscore(rast[])
        envVariableNscore = rast
        envVariableNscore[] = rastScore$nscore
        mask = envVariableNscore * 0 + 1
        mask.grid = as(mask, "SpatialGridDataFrame")
        simulation = predict(variogramModel, newdata = mask.grid, 
            nsim = 1)
        simRaster = raster(simulation) * mask
        simRaster = nscoreBack(simRaster[], rastScore)
    }
    if (quantileTransformation == TRUE) {
        mask = rast * 0 + 1
        mask.grid = as(mask, "SpatialGridDataFrame")
        simulation = predict(variogramModel, newdata = mask.grid, 
            nsim = 1)
        simRaster = raster(simulation) * mask
        envVariableData = as.data.frame(rast)
        names(envVariableData) = "value"
        envVariableData = envVariableData[with(envVariableData, 
            order(value)), ]
        simRasterData = as.data.frame(simRaster)
        names(simRasterData) = "value"
        indexes = 1:dim(simRasterData)[1]
        simRasterData$indexes = indexes
        simRasterData = simRasterData[with(simRasterData, order(value)), 
            ]
        simRasterData$value = envVariableData
        simRasterData = simRasterData[with(simRasterData, order(indexes)), 
            ]
        simRaster = simRasterData[, "value"]
    }
    rast[] = simRaster
    return(rast)
}
