spreadGraphic2 <-
function (localTreesDirectory, nberOfExtractionFiles, prob = 0.95, 
    startDatum, precision = 1) 
{
    percentage = gsub("0\\.", "", as.character(prob))
    timeInterval = precision
    nodes = c()
    for (i in 1:nberOfExtractionFiles) {
        tab = read.csv(paste0(localTreesDirectory, "/TreeExtractions_", 
            i, ".csv"), header = T)
        startingNodes = tab[, c("startYear", "startLon", "startLat")]
        colnames(startingNodes) = c("time", "lon", "lat")
        endingNodes = tab[, c("endYear", "endLon", "endLat")]
        colnames(endingNodes) = c("time", "lon", "lat")
        nodes = rbind(nodes, startingNodes, endingNodes)
        if (i == 1) 
            endDatum = max(tab[, "endYear"])
    }
    timeSlices = seq(startDatum, endDatum, timeInterval)
    spreads = list()
    c = 0
    for (i in 1:length(timeSlices)) {
        startTime = timeSlices[i]
        endTime = timeSlices[i] + timeInterval
        selectedNodes = nodes[which((nodes[, "time"] >= startTime) & 
            (nodes[, "time"] < endTime)), ]
        selectedNodes = unique(selectedNodes[, c("lon", "lat")])
        if (dim(selectedNodes)[1] > 2) {
            c = c + 1
            H = Hpi(cbind(selectedNodes[, "lon"], selectedNodes[, 
                "lat"]))
            print(timeSlices[i])
            kde = kde(cbind(selectedNodes[, "lon"], selectedNodes[, 
                "lat"]), H = H, compute.cont = T, gridsize = c(1000, 
                1000))
            contourLevel = contourLevels(kde, prob = 0.01)
            polygons = list()
            contourLines = contourLines(kde$eval.points[[1]], 
                kde$eval.points[[2]], kde$estimate, level = contourLevel)
            for (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x, 
                contourLines[[j]]$y))
            ps = Polygons(polygons, 1)
            contourPolygons = SpatialPolygons(list(ps))
            spdf = SpatialPolygonsDataFrame(contourPolygons, 
                data.frame(ID = 1:length(contourPolygons)))
            names(spdf) = round(timeSlices[i], 3)
            spreads[[c]] = spdf
        }
    }
    return(spreads)
}