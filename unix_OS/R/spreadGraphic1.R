spreadGraphic1 <-
function(localTreesDirectory, nberOfExtractionFiles, rast, prob=0.95, startDatum, precision=1, timeLayers=FALSE, nberOfCores=1, origin=FALSE, yearLayers=FALSE) {

\tshowingPlots = FALSE
\tprob = 1-prob
\tdigits = 0
\tif (precision < 1) digits = 1
\tif (precision < 0.1) digits = 2
\tif (precision < 0.01) digits = 3
\tif (precision < 0.001) digits = 4
\tif (precision < 0.0001) digits = 5
\tif (precision < 0.00001) digits = 6
\tregisterDoMC(cores=nberOfCores)
\trasts = list() # "timeLayers" option
\trast1 = rast; rast1[!is.na(rast1[])] = 1
\trast2 = rast; rast2[!is.na(rast2[])] = NA
\teA = extent(xmin(rast), xmax(rast), ymin(rast), ymax(rast))
\tfileName = paste(localTreesDirectory, "/TreeExtractions_1.csv", sep="")
\tdata = read.csv(fileName, h=T)
\tmostRecentSamplingDatum = max(data[,"endYear"])
\toriginPoints = matrix(nrow=nberOfExtractionFiles, ncol=3)
\txMin = min(min(data[,"startLon"]),min(data[,"endLon"]))
\txMax = max(max(data[,"startLon"]),max(data[,"endLon"]))
\tyMin = min(min(data[,"startLat"]),min(data[,"endLat"]))
\tyMax = max(max(data[,"startLat"]),max(data[,"endLat"]))
\tpoints = c(); c = 0
\tfor (t in 1:nberOfExtractionFiles)
\t\t{
\t\t\tif (mostRecentSamplingDatum < max(data[,"endYear"])) mostRecentSamplingDatum = max(data[,"endYear"])
\t\t\tfileName = paste(localTreesDirectory, "/TreeExtractions_", t, ".csv", sep="")
\t\t\tdata = read.csv(fileName, h=T)
\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\tif (timeLayers == TRUE) originPoints[t,1] = ceiling(data[data[,"startYear"]==min(data[,"startYear"]),"startYear"][1])
\t\t\tif (timeLayers == FALSE) originPoints[t,1] = data[data[,"startYear"]==min(data[,"startYear"]),"startYear"][1]
\t\t\toriginPoints[t,2] = data[data[,"startYear"]==min(data[,"startYear"]),"startLon"][1]
\t\t\toriginPoints[t,3] = data[data[,"startYear"]==min(data[,"startYear"]),"startLat"][1]
\t\t\tfor (i in 1:dim(data)[1])
\t\t\t\t{
\t\t\t\t\tc = c+1; points = rbind(points, cbind(data[i,"endYear"], data[i,"endLon"], data[i,"endLat"]))
\t\t\t\t}
\t\t\tif (xMin > min(data[,"startLon"])) xMin = min(data[,"startLon"])
\t\t\tif (xMin > min(data[,"endLon"])) xMin = min(data[,"endLon"])
\t\t\tif (xMax < min(data[,"startLon"])) xMax = min(data[,"startLon"])
\t\t\tif (xMax < min(data[,"endLon"])) xMax = min(data[,"endLon"])
\t\t\tif (yMin > min(data[,"startLat"])) yMin = min(data[,"startLat"])
\t\t\tif (yMin > min(data[,"endLat"])) yMin = min(data[,"endLat"])
\t\t\tif (yMax < min(data[,"startLat"])) yMax = min(data[,"startLat"])
\t\t\tif (yMax < min(data[,"endLat"])) yMax = min(data[,"endLat"])
\t\t}
\t
\tif (origin == TRUE)
\t\t{
\t\t\tH = Hpi(originPoints[,2:3])
\t\t\txMin = c(rast2@extent@xmin, rast2@extent@ymin)
\t\t\txMax = c(rast2@extent@xmax, rast2@extent@ymax)
\t\t\tkde = kde(originPoints[,2:3], H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
\t\t\tthreshold = contourLevels(kde, prob)
\t\t\ttemp = raster(kde); crs(temp) = crs(rast2)
\t\t\ttemp[temp[]<threshold] = NA
\t\t\ttemp[!is.na(temp[])] = i
\t\t\trast2 = raster::merge(rast2, temp)
\t\t}\telse\t{
\t\t\tpoints = rbind(originPoints, points)
\t\t\tpoints = points[order(points[,1]),]
\t\t\tstartYear = min(points[,1])
\t\t\tendYear = max(points[,1])
\t\t\trast1[rast1[]==1] = startYear
\t\t\tpts = points[points[,1]==startYear,2:3]
\t\t\tif (nberOfCores == 1)
\t\t\t\t{
\t\t\t\t\tc = 0
\t\t\t\t\t# for (i in seq(startDatum, endYear, by=precision))
\t\t\t\t\tfor (i in seq(floor(startDatum),floor(endYear),by=precision)) # NEW
\t\t\t\t\t\t{
\t\t\t\t\t\t\tpts = points[(points[,1]>=i)&(points[,1]<(i+precision)),2:3]
\t\t\t\t\t\t\tpts = unique(pts)
\t\t\t\t\t\t\tif (length(pts) > 0)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (length(pts) > 2)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (length(pts) > 4)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tH = Hpi(pts)
\t\t\t\t\t\t\t\t\t\t\t\t\txMin = c(rast2@extent@xmin, rast2@extent@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\txMax = c(rast2@extent@xmax, rast2@extent@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\tkde = kde(pts, H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp = raster(kde); crs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\t\t\tthreshold = contourLevels(kde, prob)
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp[temp[]<threshold] = NA
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp[!is.na(temp[])] = i
\t\t\t\t\t\t\t\t\t\t\t\t\trast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\tif (timeLayers == TRUE) rasts[[c]] = rast2
\t\t\t\t\t\t\t\t\t\t\t\t\tif (yearLayers == TRUE) rasts[[c]] = temp
\t\t\t\t\t\t\t\t\t\t\t\t\tif ((nberOfCores == 1) & (showingPlots == TRUE))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tplotRaster(rast=rast2); points(pts, pch=3, cex=1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttitle(main=i, cex.main=0.9, line=-1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprint(i)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp = rasterize(pts, rast2, i)
\t\t\t\t\t\t\t\t\t\t\t\t\tcrs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\t\t\trast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\tif (timeLayers == TRUE) rasts[[c]] = rast2
\t\t\t\t\t\t\t\t\t\t\t\t\tif (yearLayers == TRUE) rasts[[c]] = temp
\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tp = data.frame()
\t\t\t\t\t\t\t\t\t\t\tp[1,1] = pts[1]; p[1,2] = pts[2]
\t\t\t\t\t\t\t\t\t\t\ttemp = rasterize(p, rast2, i)
\t\t\t\t\t\t\t\t\t\t\tcrs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\trast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\tif (timeLayers == TRUE) rasts[[c]] = rast2
\t\t\t\t\t\t\t\t\t\t\tif (yearLayers == TRUE) rasts[[c]] = temp
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}\t\t\t
\t\t\t\t\t\t}
\t\t\t\t}\telse\t{
\t\t\t\t\tbuffer = list(); emptyList = list()
\t\t\t\t\tbuffer = foreach(i = seq(startDatum,(endYear+precision),by=precision)) %dopar% {\t\t\t
\t\t\t\t\t# for (i in 1:seq(startDatum,(endYear+precision),by=precision)) {
\t\t\t\t\t\t\tpts = points[(points[,1]>=i)&(points[,1]<(i+precision)),2:3]
\t\t\t\t\t\t\tpts = unique(pts)
\t\t\t\t\t\t\tif (length(pts) > 0)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (length(pts) > 2)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (length(pts) > 4)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tH = Hpi(pts)
\t\t\t\t\t\t\t\t\t\t\t\t\txMin = c(rast2@extent@xmin, rast2@extent@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\txMax = c(rast2@extent@xmax, rast2@extent@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\tkde = kde(pts, H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp = raster(kde); crs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\t\t\tthreshold = contourLevels(kde, prob)
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp[temp[]<threshold] = NA
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp[!is.na(temp[])] = i
\t\t\t\t\t\t\t\t\t\t\t\t\trast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t\t\t\tif ((nberOfCores == 1) & (showingPlots == TRUE))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tplotRast(rast=rast2)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprint(i)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\ttemp = rasterize(pts, rast2, i)
\t\t\t\t\t\t\t\t\t\t\t\t\tcrs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\t\t\t# rast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\tp = data.frame()
\t\t\t\t\t\t\t\t\t\t\tp[1,1] = pts[1]; p[1,2] = pts[2]
\t\t\t\t\t\t\t\t\t\t\ttemp = rasterize(p, rast2, i)
\t\t\t\t\t\t\t\t\t\t\tcrs(temp) = crs(rast2)
\t\t\t\t\t\t\t\t\t\t\t# rast2 = raster::merge(rast2, temp)
\t\t\t\t\t\t\t\t\t\t}\t\t\t
\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\ttemp = "empty" 
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t# buffer[[i]] = temp
\t\t\t\t\t\t\ttemp
\t\t\t\t\t\t}
\t\t\t\t\tc = 0
\t\t\t\t\tfor (i in 1:length(buffer))
\t\t\t\t\t\t{\t
\t\t\t\t\t\t\tif (!is.character(buffer[[i]]))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\trast2 = raster::merge(rast2, buffer[[i]])
\t\t\t\t\t\t\t\t\tif (timeLayers == TRUE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\trasts[[c]] = rast2
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tif (yearLayers == TRUE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\trasts[[c]] = buffer[[i]]
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}\t\t\t
\t\t\tif (showingPlots == TRUE)
\t\t\t\t{
\t\t\t\t\trast1[rast1[]==1] = min(rast2[]); rast2 = raster::merge(rast2, rast1); rast2[is.na(rast1[])] = NA
\t\t\t\t\tcols = colorRampPalette(c("grey90","red"),bias=1)(endYear-startYear)[1:(endYear-startYear)]
\t\t\t\t\tplotRast(rast=rast2, cols=cols)
\t\t\t\t}
\t\t}
\tif ((timeLayers == FALSE)&(yearLayers == FALSE))
\t\t{
\t\t\treturn(rast2)
\t\t}\telse\t{
\t\t\treturn(rasts)
\t\t}
}
