spreadGraphic1 = function(localTreesDirectory, nberOfExtractionFiles, rast, prob=0.95, startDatum, precision=1, 
						  timeLayers=FALSE, nberOfCores=1, origin=FALSE, yearLayers=FALSE) {

	showingPlots = FALSE
	prob = 1-prob
	digits = 0
	if (precision < 1) digits = 1
	if (precision < 0.1) digits = 2
	if (precision < 0.01) digits = 3
	if (precision < 0.001) digits = 4
	if (precision < 0.0001) digits = 5
	if (precision < 0.00001) digits = 6
	# registerDoMC(cores=nberOfCores)
	rasts = list() # "timeLayers" option
	rast1 = rast; rast1[!is.na(rast1[])] = 1
	rast2 = rast; rast2[!is.na(rast2[])] = NA
	eA = extent(xmin(rast), xmax(rast), ymin(rast), ymax(rast))
	fileName = paste(localTreesDirectory,"/TreeExtractions_1.csv",sep="")
	data = read.csv(fileName, h=T)
	mostRecentSamplingDatum = max(data[,"endYear"])
	originPoints = matrix(nrow=nberOfExtractionFiles, ncol=3)
	xMin = min(min(data[,"startLon"]),min(data[,"endLon"]))
	xMax = max(max(data[,"startLon"]),max(data[,"endLon"]))
	yMin = min(min(data[,"startLat"]),min(data[,"endLat"]))
	yMax = max(max(data[,"startLat"]),max(data[,"endLat"]))
	points = c(); c = 0
	for (t in 1:nberOfExtractionFiles)
		{
			if (mostRecentSamplingDatum < max(data[,"endYear"])) mostRecentSamplingDatum = max(data[,"endYear"])
			fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear, startYear)),]
			if (timeLayers == TRUE) originPoints[t,1] = ceiling(data[data[,"startYear"]==min(data[,"startYear"]),"startYear"][1])
			if (timeLayers == FALSE) originPoints[t,1] = data[data[,"startYear"]==min(data[,"startYear"]),"startYear"][1]
			originPoints[t,2] = data[data[,"startYear"]==min(data[,"startYear"]),"startLon"][1]
			originPoints[t,3] = data[data[,"startYear"]==min(data[,"startYear"]),"startLat"][1]
			for (i in 1:dim(data)[1])
				{
					c = c+1; points = rbind(points, cbind(data[i,"endYear"], data[i,"endLon"], data[i,"endLat"]))
				}
			if (xMin > min(data[,"startLon"])) xMin = min(data[,"startLon"])
			if (xMin > min(data[,"endLon"])) xMin = min(data[,"endLon"])
			if (xMax < min(data[,"startLon"])) xMax = min(data[,"startLon"])
			if (xMax < min(data[,"endLon"])) xMax = min(data[,"endLon"])
			if (yMin > min(data[,"startLat"])) yMin = min(data[,"startLat"])
			if (yMin > min(data[,"endLat"])) yMin = min(data[,"endLat"])
			if (yMax < min(data[,"startLat"])) yMax = min(data[,"startLat"])
			if (yMax < min(data[,"endLat"])) yMax = min(data[,"endLat"])
		}
	
	if (origin == TRUE)
		{
			H = Hpi(originPoints[,2:3])
			xMin = c(rast2@extent@xmin, rast2@extent@ymin)
			xMax = c(rast2@extent@xmax, rast2@extent@ymax)
			kde = kde(originPoints[,2:3], H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
			threshold = contourLevels(kde, prob)
			temp = raster(kde); crs(temp) = crs(rast2)
			temp[temp[]<threshold] = NA
			temp[!is.na(temp[])] = i
			rast2 = raster::merge(rast2, temp)
		}	else	{
			points = rbind(originPoints, points)
			points = points[order(points[,1]),]
			startYear = min(points[,1])
			endYear = max(points[,1])
			rast1[rast1[]==1] = startYear
			pts = points[points[,1]==startYear,2:3]
			if (nberOfCores == 1)
				{
					c = 0
					# for (i in seq(startDatum, endYear, by=precision))
					for (i in seq(floor(startDatum),floor(endYear),by=precision)) # NEW
						{
							pts = points[(points[,1]>=i)&(points[,1]<(i+precision)),2:3]
							pts = unique(pts)
							if (length(pts) > 0)
								{
									if (length(pts) > 2)
										{
											if (length(pts) > 4)
												{
													H = Hpi(pts)
													xMin = c(rast2@extent@xmin, rast2@extent@ymin)
													xMax = c(rast2@extent@xmax, rast2@extent@ymax)
													kde = kde(pts, H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
													temp = raster(kde); crs(temp) = crs(rast2)
													threshold = contourLevels(kde, prob)
													temp[temp[]<threshold] = NA
													temp[!is.na(temp[])] = i
													rast2 = raster::merge(rast2, temp)
													c = c+1
													if (timeLayers == TRUE) rasts[[c]] = rast2
													if (yearLayers == TRUE) rasts[[c]] = temp
													if ((nberOfCores == 1) & (showingPlots == TRUE))
														{
															plotRaster(rast=rast2); points(pts, pch=3, cex=1)
															title(main=i, cex.main=0.9, line=-1)
														}	else	{
															print(i)
														}
												}	else	{
													temp = rasterize(pts, rast2, i)
													crs(temp) = crs(rast2)
													rast2 = raster::merge(rast2, temp)
													c = c+1
													if (timeLayers == TRUE) rasts[[c]] = rast2
													if (yearLayers == TRUE) rasts[[c]] = temp
												}	
										}	else	{
											p = data.frame()
											p[1,1] = pts[1]; p[1,2] = pts[2]
											temp = rasterize(p, rast2, i)
											crs(temp) = crs(rast2)
											rast2 = raster::merge(rast2, temp)
											c = c+1
											if (timeLayers == TRUE) rasts[[c]] = rast2
											if (yearLayers == TRUE) rasts[[c]] = temp
										}
								}			
						}
				}	else	{
					buffer = list(); emptyList = list()
					# buffer = foreach(i = seq(startDatum,(endYear+precision),by=precision)) %dopar% {			
					for (i in 1:seq(startDatum,(endYear+precision),by=precision)) {
							pts = points[(points[,1]>=i)&(points[,1]<(i+precision)),2:3]
							pts = unique(pts)
							if (length(pts) > 0)
								{
									if (length(pts) > 2)
										{
											if (length(pts) > 4)
												{
													H = Hpi(pts)
													xMin = c(rast2@extent@xmin, rast2@extent@ymin)
													xMax = c(rast2@extent@xmax, rast2@extent@ymax)
													kde = kde(pts, H=H, gridsize=c(rast2@ncols, rast2@nrows), xmin=xMin, xmax=xMax)
													temp = raster(kde); crs(temp) = crs(rast2)
													threshold = contourLevels(kde, prob)
													temp[temp[]<threshold] = NA
													temp[!is.na(temp[])] = i
													rast2 = raster::merge(rast2, temp)
													if ((nberOfCores == 1) & (showingPlots == TRUE))
														{
															plotRast(rast=rast2)
														}	else	{
															print(i)
														}
												}	else	{
													temp = rasterize(pts, rast2, i)
													crs(temp) = crs(rast2)
													# rast2 = raster::merge(rast2, temp)
												}	
										}	else	{
											p = data.frame()
											p[1,1] = pts[1]; p[1,2] = pts[2]
											temp = rasterize(p, rast2, i)
											crs(temp) = crs(rast2)
											# rast2 = raster::merge(rast2, temp)
										}			
								}	else	{
									temp = "empty" 
								}
							buffer[[i]] = temp
							# temp
						}
					c = 0
					for (i in 1:length(buffer))
						{	
							if (!is.character(buffer[[i]]))
								{
									rast2 = raster::merge(rast2, buffer[[i]])
									if (timeLayers == TRUE)
										{
											c = c+1
											rasts[[c]] = rast2
										}
									if (yearLayers == TRUE)
										{
											c = c+1
											rasts[[c]] = buffer[[i]]
										}
								}
						}
				}			
			if (showingPlots == TRUE)
				{
					rast1[rast1[]==1] = min(rast2[]); rast2 = raster::merge(rast2, rast1); rast2[is.na(rast1[])] = NA
					cols = colorRampPalette(c("grey90","red"),bias=1)(endYear-startYear)[1:(endYear-startYear)]
					plotRast(rast=rast2, cols=cols)
				}
		}
	if ((timeLayers == FALSE)&(yearLayers == FALSE))
		{
			return(rast2)
		}	else	{
			return(rasts)
		}
}
