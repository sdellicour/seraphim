treesRandomisation = function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), randomProcedure=3, showingPlots=F) {

	nberOfRandomisations = 1
	# registerDoMC(cores=nberOfCores)
	hull_polygons = list()
	rotation = function(pt1, pt2, angle)
		{
			s = sin(angle); c = cos(angle)
			x = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
			x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
			x_new = x_new+pt1[1]; y_new = y_new+pt1[2]
			return(c(x_new,y_new))
		}
	nullRaster = envVariables[[1]]
	nullRaster[!is.na(nullRaster)] = 1
	names(nullRaster) = "null_raster"
	envVariables0 = envVariables	
	newEnvVariables = list(nullRaster)
	for (h in 1:length(envVariables0))
		{
			newEnvVariables[[h+1]] = envVariables0[[h]]
			newEnvVariables[[h+1]][newEnvVariables[[h+1]]<0] = NA
		}
	envVariables = newEnvVariables
	branchRandomisation3 = FALSE
	branchRandomisation2 = FALSE
	branchRandomisation1 = FALSE
	if (randomProcedure == 3) branchRandomisation3 = TRUE
	if (randomProcedure == 4)
		{
			branchRandomisation2 = TRUE; rotatingEndNodes = TRUE
		}
	if (randomProcedure == 5)
		{
			branchRandomisation2 = TRUE; rotatingEndNodes = FALSE
		}
	if (randomProcedure == 6) branchRandomisation1 = TRUE
	nberOfConnections = matrix(nrow=1, ncol=nberOfExtractionFiles)	
	totalnberOfConnections = 0
	extractionFileName = "TreeExtractions"
	if (nchar(localTreesDirectory) == 0)
		{
			data = read.csv(paste(extractionFileName,"_1.csv",sep=""), header=T, dec=".")
		}	else	{
			data = read.csv(paste(localTreesDirectory,"/",extractionFileName,"_1.csv",sep=""), header=T, dec=".")	
		}
	data = data[with(data, order(startYear,endYear)),]
	node1 = list()
	node2 = list()
	startYear = list()
	dispersalTime = list()
	treeIDs = list()
	dispersalRate = list()
	fromCoor = list()
	toCoor = list()
	datas = list()
	for (t in 1:nberOfExtractionFiles)
		{
			if (t != 1)
				{
					if (nchar(localTreesDirectory) == 0)
						{
							fileName = paste(extractionFileName,"_",t,".csv",sep="")
						}	else	{
							fileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
						}	
					data = read.csv(fileName, h = T)
					data = data[with(data, order(endYear, startYear)),]
				}
			ancestralNodeNAonNullRaster = TRUE
			while (ancestralNodeNAonNullRaster == TRUE)
				{
					ancestralNodeNAonNullRaster = FALSE
					ancestralBranches = which(!data[,"node1"]%in%data[,"node2"])
					indicesOfBranchesToRemove = c()
					for (i in 1:length(ancestralBranches))
						{
							if (is.na(raster::extract(nullRaster, cbind(data[ancestralBranches[i],"startLon"],data[ancestralBranches[i],"startLat"]))))
								{
									ancestralNodeNAonNullRaster = TRUE
									indicesOfBranchesToRemove = c(indicesOfBranchesToRemove, ancestralBranches[i])
								}
						}
					if (length(indicesOfBranchesToRemove) > 0)
						{
							data = data[-indicesOfBranchesToRemove,]
						}
				}
			nberOfConnections[t] = dim(data)[1]		
			node1[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)	
			node1[[t]][] = data[,"node1"]
			node2[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)	
			node2[[t]][] = data[,"node2"]
			startYear[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)	
			startYear[[t]][] = data[,"startYear"]
			dispersalTime[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)
			dispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
			colnames(dispersalTime[[t]]) = "dispersalTime"
			fromCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
			fromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
			toCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
			toCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
			totalnberOfConnections = totalnberOfConnections + nberOfConnections[t] 
			if (("treeID"%in%colnames(data)) == TRUE)
				{
					treeIDs[[t]] = data[1,"treeID"]
				}	else	{
					treeIDs[[t]] = "noTreeID"
				}
			datas[[t]] = data
		}
	hullRasters = list()
	hullRasters[1:length(envVariables)] = envVariables[1:length(envVariables)]
	points = matrix(nrow=(totalnberOfConnections*2),ncol=2)
	a = 0
	for (t in 1:nberOfExtractionFiles)
		{
			if (t > 1)
				{
					a = a + nberOfConnections[t-1]
				}
			for (i in 1:nberOfConnections[t])
				{
					index = (a*2) + ((i-1)*2) + 1
					points[index,1] = fromCoor[[t]][i,1]
					points[index,2] = fromCoor[[t]][i,2]
					points[(index+1),1] = toCoor[[t]][i,1]
					points[(index+1),2] = toCoor[[t]][i,2]								
				}
		}
	points = points[points[,1]>extent(hullRasters[[1]])@xmin,]
	points = points[points[,1]<extent(hullRasters[[1]])@xmax,]
	points = points[points[,2]>extent(hullRasters[[1]])@ymin,]
	points = points[points[,2]<extent(hullRasters[[1]])@ymax,]
	if (length(hull_polygons) == 0)
		{		
			hull = chull(points)
			hull = c(hull,hull[1])
			# plot(points); lines(points[hull,])
			p = Polygon(points[hull,])
			ps = Polygons(list(p),1)
			sps = SpatialPolygons(list(ps))
			# c = mask(envVariables[[2]],sps); plot(c)
		}
	if (length(hull_polygons) > 0)
		{
			sps = hull_polygons
		}
	pointsRaster = rasterize(points, crop(hullRasters[[1]], sps, snap="out"))
	pointsRaster[!is.na(pointsRaster[])] = 0		
	for (h in 1:length(envVariables))
		{
			# plot(mask(simRasters[[h]],sps))
			hullRasters[[h]] = crop(hullRasters[[h]], sps, snap="out")
			bufferRaster = hullRasters[[h]]
			hullRasters[[h]] = mask(hullRasters[[h]], sps, snap="out")
			hullRasters[[h]][!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
			names(hullRasters[[h]]) = gsub(".asc","",names(envVariables[[h]]))
			names(hullRasters[[h]]) = gsub(".tif","",names(envVariables[[h]]))
			names(hullRasters[[h]]) = gsub(".gri","",names(envVariables[[h]]))
		}
	if (nberOfRandomisations > 0)
		{	
			for (s in 1:nberOfRandomisations)
				{
					if (branchRandomisation3 == TRUE)
						{
							simRasters = list(); simRasters = hullRasters
							fromCoorRand = fromCoor; toCoorRand = toCoor
							buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									if (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
										{
											counter1 = 0; # print(t)
											twoPointsOnTheGrid = FALSE
											while (twoPointsOnTheGrid == FALSE)
												{
													twoPointsOnTheGrid = TRUE; counter1 = counter1+1
													if (counter1 == 1)
														{
															cat("Randomising tree ",t,"\n",sep="")
														}	else		{
															cat("Randomising tree ",t,", again","\n",sep="")
														}
													fromCoorRand[[t]][,] = NA; toCoorRand[[t]][,] = NA
													if (showingPlots == TRUE)
														{
															if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
															if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
															lines(points[hull,], lwd=0.5, col="black")
															text1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
															mtext(text1, col="black", cex=0.7, line=0)
															# for (i in 1:dim(fromCoor[[t]])[1])
																# {
																	# segments(fromCoor[[t]][i,1], fromCoor[[t]][i,2], toCoor[[t]][i,1], toCoor[[t]][i,2], col="black", lwd=0.3)
																# }
														}
													ancestralIndex = list(); ancestralNodes = list(); counter = 0
													for (i in 1:length(node1[[t]]))
														{
															ancestralNodeBoolean = TRUE
															for (j in 1:length(node2[[t]]))
																{
																	if (node1[[t]][i,1] == node2[[t]][j,1])
																		{
																			ancestralNodeBoolean = FALSE
																		}
																}
															if (ancestralNodeBoolean == TRUE)
																{
																	counter = counter+1
																	ancestralIndex[[counter]] = i
																	ancestralNodes[[counter]] = node1[[t]][i,1]
																}	
														}
													for (i in 1:length(ancestralIndex))
														{
															fromCoorRand[[t]][ancestralIndex[[i]],1] = fromCoor[[t]][ancestralIndex[[i]],1]
															fromCoorRand[[t]][ancestralIndex[[i]],2] = fromCoor[[t]][ancestralIndex[[i]],2]
														}
													ancestralNodes = unique(ancestralNodes)
													startingNodes = list(); startingNodes = ancestralNodes # startingNodes[[1]] = ancestralNode[1]
													while (length(startingNodes) > 0)
														{
															newStartingNodes = list(); c = 0
															for (i in 1:length(startingNodes))
																{	
																	nodes2 = node2[[t]][which(node1[[t]][,1]==startingNodes[[i]]),1]
																	if (length(nodes2) > 0)
																		{
																			for (j in 1:length(nodes2))
																				{
																					c = c+1
																					newStartingNodes[[c]] = nodes2[j]				
																					k = which(node2[[t]][,1]==nodes2[j])
																					pt01 = c(fromCoor[[t]][k,1], fromCoor[[t]][k,2])
																					pt02 = c(toCoor[[t]][k,1], toCoor[[t]][k,2])
																					xTranslation = pt02[1]-pt01[1]
																					yTranslation = pt02[2]-pt01[2]
																					pt1 = c(fromCoorRand[[t]][k,1], fromCoorRand[[t]][k,2])
																					pt2 = c(NA,NA)
																					pt2[1] = pt1[1]+xTranslation
																					pt2[2] = pt1[2]+yTranslation
																					pt2NAarea = TRUE
																					counter2 = 0
																					while (pt2NAarea == TRUE)
																						{
																							counter2 = counter2+1	
																							onTheGrid = TRUE
																							angle = (2*pi)*runif(1)
																							pt2_rotated = rotation(pt1, pt2, angle)
																							if (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
																								{
																									onTheGrid = FALSE
																								}
																							if (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
																								{
																									onTheGrid = FALSE
																								}	
																							if (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
																								{
																									onTheGrid = FALSE
																								}
																							if (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
																								{
																									onTheGrid = FALSE
																								}
																							if (onTheGrid == TRUE)
																								{
																									NAarea = FALSE
																									for (h in 1:length(hullRasters))
																										{
																											if (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
																												{
																													NAarea = TRUE
																												}
																										}
																									if (NAarea == FALSE)
																										{	
																											pt2NAarea = FALSE
																										}		
																								}									
																							if (counter2 > 100)
																								{
																									# print("counter2 > 100")
																									ancestralNodeID = FALSE
																									for (h in 1:length(ancestralNodes))
																										{
																											if (ancestralNodes[[h]] == node1[[t]][which(node2[[t]][,1]==nodes2[j]),1])
																												{
																													ancestralNodeID = TRUE
																												}
																										}
																									if (ancestralNodeID == FALSE)
																										{
																											if (counter1 <= 10)
																												{
																													# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
																													pt2_rotated = pt1
																													pt2NAarea = FALSE
																													twoPointsOnTheGrid = FALSE
																												}	else	{
																													# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
																													pt2_rotated = pt1
																													pt2NAarea = FALSE
																													twoPointsOnTheGrid = TRUE
																												}
																										}
																									if (ancestralNodeID == TRUE)
																										{
																											# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
																											pt2_rotated = pt1
																											pt2NAarea = FALSE
																											twoPointsOnTheGrid = TRUE
																										}
																								}							
																						}
																					pt2 = pt2_rotated
																					if (showingPlots == TRUE)
																						{
																							points(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
																							points(pt2[1], pt2[2], pch=16, col="black", cex=0.25)
																							segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
																						}
																					fromCoorRand[[t]][k,1] = pt1[1]; fromCoorRand[[t]][k,2] = pt1[2]
																					toCoorRand[[t]][k,1] = pt2[1]; toCoorRand[[t]][k,2] = pt2[2]
																					toModify = which(node1[[t]][,1]==nodes2[j])
																					for (k in 1:length(toModify))
																						{
																							fromCoorRand[[t]][toModify[k],1] = pt2[1]
																							fromCoorRand[[t]][toModify[k],2] = pt2[2]
																						}
																				}
																		}
																}		
															startingNodes = newStartingNodes	
														}
												}
											temp = datas[[t]]
											temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
											write.csv(temp, paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
										}
									# t
								}
						}
					if (branchRandomisation2 == TRUE)
						{
							simRasters = list(); simRasters = hullRasters
							cat("Analysis of randomised branch positions ", s, "\n", sep="")
							fromCoorRand = fromCoor; toCoorRand = toCoor
							buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									if (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
										{
											if (showingPlots == TRUE)
												{
													if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
													if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
													lines(points[hull,], lwd=0.5, col="black")
													text1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
													mtext(text1, col="black", cex=0.7, line=0)
												}
											for (i in 1:nberOfConnections[t])
												{
													if (rotatingEndNodes == TRUE)
														{	
															pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
															pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
														}	else	{
															pt1 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
															pt2 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
														}
													twoPointsOnTheGrid = FALSE
													while (twoPointsOnTheGrid == FALSE)
														{	
															pt2NAarea = TRUE
															counter = 0
															while (pt2NAarea == TRUE)
																{
																	counter = counter+1	
																	onTheGrid = TRUE
																	angle = (2*pi)*runif(1)
																	pt2_rotated = rotation(pt1, pt2, angle)
																	if (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
																		{
																			onTheGrid = FALSE
																		}
																	if (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
																		{
																			onTheGrid = FALSE
																		}	
																	if (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
																		{
																			onTheGrid = FALSE
																		}
																	if (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
																		{
																			onTheGrid = FALSE
																		}
																	if (onTheGrid == TRUE)
																		{
																			NAarea = FALSE
																			for (h in 1:length(hullRasters))
																				{
																					if (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
																						{
																							NAarea = TRUE
																							twoPointsOnTheGrid = TRUE
																						}
																				}
																			if (NAarea == FALSE)
																				{	
																					pt2NAarea = FALSE
																					twoPointsOnTheGrid = TRUE
																				}		
																		}
																	if (counter > 100)
																		{
																			pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																			pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
																			pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
																		}							
																}
														}
													pt2 = pt2_rotated
													if (showingPlots == TRUE)
														{
															points(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
															points(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
															segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
														}
													if (rotatingEndNodes == TRUE)
														{
															fromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
															toCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
														}	else	{
															toCoorRand[[t]][i,1] = pt1[1]; toCoorRand[[t]][i,2] = pt1[2]
															fromCoorRand[[t]][i,1] = pt2[1]; fromCoorRand[[t]][i,2] = pt2[2]
														}
												}
											temp = datas[[t]]
											temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
											write.csv(temp, paste0(localTreeDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
										}
									# t
								}
						}			
					if (branchRandomisation1 == TRUE)
						{
							simRasters = list(); simRasters = hullRasters
							cat("Analysis of randomised branch positions ", s, "\n", sep="")
							fromCoorRand = fromCoor; toCoorRand = toCoor
							buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									if (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
										{
											if (showingPlots == TRUE)
												{
													if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
													if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=T)
													plot(sps, lwd=0.5, border="black", add=T)
													text1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
													mtext(text1, col="black", cex=0.7, line=0)
												}
											for (i in 1:nberOfConnections[t])
												{
													pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
													pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
													twoPointsOnTheGrid = FALSE
													while (twoPointsOnTheGrid == FALSE)
														{	
															pt1NAarea = TRUE
															while (pt1NAarea == TRUE)
																{
																	pt1_translated = pt1
																	xTranslation = runif(1)*(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
																	yTranslation = runif(1)*(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
																	pt1_translated[1] = pt1[1]+xTranslation; pt1_translated[2] = pt1[2]+yTranslation
																	if (pt1_translated[1] > extent(hullRasters[[1]])@xmax)
																		{
																			pt1_translated[1] = pt1_translated[1]-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
																			xTranslation = xTranslation-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
																		}
																	if (pt1_translated[2] > extent(hullRasters[[1]])@ymax)
																		{
																			pt1_translated[2] = pt1_translated[2]-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
																			yTranslation = yTranslation-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
																		}
																	NAarea = FALSE	
																	for (h in 1:length(hullRasters))
																		{
																			if (is.na(raster::extract(hullRasters[[h]],cbind(pt1_translated[1],pt1_translated[2]))))
																				{
																					NAarea = TRUE
																				}
																		}
																	insideAtLeastOneHullPolygon = FALSE
																	for (h in 1:length(sps))
																		{
																			pol_x = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,1]
																			pol_y = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,2]
																			if (point.in.polygon(pt1_translated[1], pt1_translated[2], pol_x, pol_y) == 1)
																				{
																					insideAtLeastOneHullPolygon = TRUE
																				}
																		}
																	if (insideAtLeastOneHullPolygon == FALSE) NAarea = TRUE
																	if (NAarea == FALSE)
																		{
																			pt1NAarea = FALSE
																			pt1 = pt1_translated
																			pt2[1] = pt2[1]+xTranslation; pt2[2] = pt2[2]+yTranslation
																		}			
																}
															pol_index = NA
															for (j in 1:length(sps))
																{
																	pol_x = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,1]
																	pol_y = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,2]
																	if (point.in.polygon(pt1[1], pt1[2], pol_x, pol_y) == 1)
																		{
																			pol_index = j
																		}
																}
															pt2NAarea = TRUE
															counter = 0
															while (pt2NAarea == TRUE)
																{
																	counter = counter+1	
																	onTheGrid = TRUE
																	angle = (2*pi)*runif(1)
																	pt2_rotated = rotation(pt1, pt2, angle)
																	if (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
																		{
																			onTheGrid = FALSE
																		}
																	if (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
																		{
																			onTheGrid = FALSE
																		}	
																	if (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
																		{
																			onTheGrid = FALSE
																		}
																	if (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
																		{
																			onTheGrid = FALSE
																		}
																	if (onTheGrid == TRUE)
																		{
																			NAarea = FALSE
																			for (h in 1:length(hullRasters))
																				{
																					if (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
																						{
																							NAarea = TRUE
																							twoPointsOnTheGrid = TRUE
																						}
																				}
																			if (NAarea == FALSE)
																				{
																					pol_x = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,1]
																					pol_y = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,2]
																					if (point.in.polygon(pt2_rotated[1], pt2_rotated[2], pol_x, pol_y) == 1)
																						{
																							pt2NAarea = FALSE
																							twoPointsOnTheGrid = TRUE
																						}
																				}		
																		}
																	if (counter > 100)
																		{
																			pt2NAarea = FALSE
																			pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
																			pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
																			twoPointsOnTheGrid = FALSE
																		}							
																}
														}
													pt2 = pt2_rotated
													if (showingPlots == TRUE)
														{
															points(pt1_translated[1], pt1_translated[2], pch=16, col="black", cex=0.25)
															points(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
															segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
														}
													fromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
													toCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
												}
											temp = datas[[t]]
											temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
											write.csv(temp, paste0(localTreeDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
										}
									# t
								}
						}	
				}
		}		
}
