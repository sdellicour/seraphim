treesRandomisation = function(localTreesDirectory="", randomisationDirectory="", nberOfExtractionFiles=1, envVariable, randomProcedure=3, 
							  resistance=NULL, overwrite=FALSE, showingPlots=FALSE) {

	if ((!is.null(repulsion))&(!is.null(resistance)))
		{
			stop("Both \"repulsion\" and \"resistance\" parameters are specified as non null - pick one of the two")
		}
	dir.create(file.path(randomisationDirectory), showWarnings=F)
	nberOfRandomisations = 1; # registerDoMC(cores=nberOfCores)
	rotation = function(pt1, pt2, angle)
		{
			s = sin(angle); c = cos(angle)
			x = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
			x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
			x_new = x_new+pt1[1]; y_new = y_new+pt1[2]
			return(c(x_new,y_new))
		}
	branchRandomisation3 = FALSE; branchRandomisation2 = FALSE; branchRandomisation1 = FALSE
		# Note: the repulsion/attraction is not implemented for the "randomProcedure = 1"
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
	nullRaster = envVariable; nullRaster[!is.na(nullRaster[])] = 1
	nberOfConnections = rep(NA, nberOfExtractionFiles); totalNberOfConnections = 0
	node1 = list(); node2 = list(); fromCoor = list(); toCoor = list(); datas = list()
	for (t in 1:nberOfExtractionFiles)
		{
			if (nchar(localTreesDirectory) == 0)
				{
					fileName = paste("TreeExtractions_",t,".csv",sep="")
				}	else	{
					fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
				}	
			data = read.csv(fileName, head=T)
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
			datas[[t]] = data; nberOfConnections[t] = dim(data)[1]		
			node1[[t]] = matrix(nrow=nberOfConnections[t], ncol=1); node1[[t]][] = data[,"node1"]
			node2[[t]] = matrix(nrow=nberOfConnections[t], ncol=1); node2[[t]][] = data[,"node2"]
			fromCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
			fromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
			toCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
			toCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
			totalNberOfConnections = totalNberOfConnections + nberOfConnections[t]
		}
	hullRaster = envVariable; points = matrix(nrow=(totalNberOfConnections*2),ncol=2); a = 0
	for (t in 1:nberOfExtractionFiles)
		{
			if (t > 1) a = a + nberOfConnections[t-1]
			for (i in 1:nberOfConnections[t])
				{
					index = (a*2) + ((i-1)*2) + 1
					points[index,1] = fromCoor[[t]][i,1]
					points[index,2] = fromCoor[[t]][i,2]
					points[(index+1),1] = toCoor[[t]][i,1]
					points[(index+1),2] = toCoor[[t]][i,2]
				}
		}
	points = points[points[,1]>extent(hullRaster)@xmin,]
	points = points[points[,1]<extent(hullRaster)@xmax,]
	points = points[points[,2]>extent(hullRaster)@ymin,]
	points = points[points[,2]<extent(hullRaster)@ymax,]
	hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
	ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
	pointsRaster = rasterize(points, crop(hullRaster, sps, snap="out"))
	pointsRaster[!is.na(pointsRaster[])] = 0		
	hullRaster = crop(hullRaster, sps, snap="out"); bufferRaster = hullRaster
	hullRaster = mask(hullRaster, sps, snap="out")
	hullRaster[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
	names(hullRaster) = gsub(".asc","",names(envVariable))
	names(hullRaster) = gsub(".tif","",names(envVariable))
	names(hullRaster) = gsub(".gri","",names(envVariable))
	if (nberOfRandomisations > 0)
		{	
			for (s in 1:nberOfRandomisations)
				{
					if (branchRandomisation3 == TRUE)
						{
							fromCoorRand = fromCoor; toCoorRand = toCoor; buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									newRandomisation = TRUE
									if (file.exists(paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv")))
										{
											newRandomisation = FALSE
											if (overwrite) newRandomisation = TRUE
										}
									if (newRandomisation)
										{
											counter1 = 0; # print(t)
											twoPointsOnTheGrid = FALSE
											while (twoPointsOnTheGrid == FALSE)
												{
													twoPointsOnTheGrid = TRUE; counter1 = counter1+1
													if (counter1 == 1)
														{
															cat("Randomising tree ",t,"\n",sep="")
														}	else	{
															cat("Randomising tree ",t,", again","\n",sep="")
														}
													fromCoorRand[[t]][,] = NA; toCoorRand[[t]][,] = NA
													if (showingPlots == TRUE)
														{
															if (t == 1) plotRaster(envVariable, addLegend=T, new=T)
															if (t >= 2) plotRaster(envVariable, addLegend=T, new=F)
															lines(points[hull,], lwd=1.5, col=rgb(222,67,39,220,maxColorValue=255)) # red
															text1 = paste("Randomisation of branch positions, sampled tree ",t,sep="")
															mtext(text1, col="gray30", cex=0.9, line=-1)
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
													startingNodes = list(); startingNodes = ancestralNodes
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
																					if (is.null(repulsion)) nberOfRotations = 1
																					if (!is.null(repulsion)) nberOfRotations = 100
																					pt2_rotated_list = list(); envValues1 = list()
																					for (g in 1:nberOfRotations)
																						{
																							pt2NAarea = TRUE; counter2 = 0
																							while (pt2NAarea == TRUE)
																								{
																									counter2 = counter2+1; onTheGrid = TRUE
																									angle = (2*pi)*runif(1); pt2_rotated = rotation(pt1, pt2, angle)
																									if (pt2_rotated[1] > extent(hullRaster)@xmax) onTheGrid = FALSE
																									if (pt2_rotated[1] < extent(hullRaster)@xmin) onTheGrid = FALSE
																									if (pt2_rotated[2] > extent(hullRaster)@ymax) onTheGrid = FALSE
																									if (pt2_rotated[2] < extent(hullRaster)@ymin) onTheGrid = FALSE
																									if (onTheGrid == TRUE)
																										{
																											NAarea = FALSE
																											if (is.na(raster::extract(hullRaster,cbind(pt2_rotated[1],pt2_rotated[2]))))
																												{
																													NAarea = TRUE
																												}
																											if (NAarea == FALSE) pt2NAarea = FALSE
																										}
																									if (counter2 > 100)
																										{
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
																															pt2_rotated = pt1; pt2NAarea = FALSE; twoPointsOnTheGrid = FALSE
																														}	else	{
																															pt2_rotated = pt1; pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																														}
																												}
																											if (ancestralNodeID == TRUE)
																												{
																													pt2_rotated = pt1; pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																												}
																										}
																								}
																							pt2_rotated_list[[g]] = pt2_rotated
																							envValues1[[g]] = raster::extract(envVariable, cbind(pt2_rotated[1],pt2_rotated[2]))
																						}
																					if (is.null(repulsion))
																						{
																							pt2 = pt2_rotated_list[[1]]
																						}	else	{
																							envValues2 = unlist(envValues1)
																							if (repulsion == TRUE)
																								{
																									buffer = envValues2-min(envValues2,na.rm=T)
																									envValues2 = max(envValues2,na.rm=T)-buffer
																								}
																							vS = envValues2-min(envValues2,na.rm=T)
																							probas = vS/sum(vS,na.rm=T)
																							probas[which(is.na(probas))] = 0
																							if (sum(probas) != 0)
																								{
																									index = sample(1:length(envValues2), 1, prob=probas)
																									pt2 = pt2_rotated_list[[index]]
																								}
																						}
																					if (showingPlots == TRUE)
																						{
																							segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.4)
																							points(pt2[1], pt2[2], pch=16, col="black", cex=0.5)
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
											write.csv(temp, paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
											if (showingPlots == TRUE) dev.copy2pdf(file=paste0(randomisationDirectory,"/TreeRandomisation_",t,".pdf"))
										}
									t
								}
						}
					if (branchRandomisation2 == TRUE)
						{
							fromCoorRand = fromCoor; toCoorRand = toCoor; buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									newRandomisation = TRUE
									if (file.exists(paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv")))
										{
											newRandomisation = FALSE
											if (overwrite) newRandomisation = TRUE
										}
									if (newRandomisation)
										{
											cat("Randomising tree ",t,"\n",sep="")
											if (showingPlots == TRUE)
												{
													if (t == 1) plotRaster(envVariable, addLegend=T, new=T)
													if (t >= 2) plotRaster(envVariable, addLegend=T, new=F)
													lines(points[hull,], lwd=1.5, col=rgb(222,67,39,220,maxColorValue=255)) # red
													text1 = paste("Randomisation of branch positions, sampled tree ",t,sep="")
													mtext(text1, col="gray30", cex=0.9, line=-1)
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
													if (is.null(repulsion)) nberOfRotations = 1
													if (!is.null(repulsion)) nberOfRotations = 100
													pt2_rotated_list = list(); envValues1 = list()
													for (g in 1:nberOfRotations)
														{
															twoPointsOnTheGrid = FALSE
															while (twoPointsOnTheGrid == FALSE)
																{	
																	pt2NAarea = TRUE; counter = 0
																	while (pt2NAarea == TRUE)
																		{
																			counter = counter+1; onTheGrid = TRUE
																			angle = (2*pi)*runif(1); pt2_rotated = rotation(pt1, pt2, angle)
																			if (pt2_rotated[1] > extent(hullRaster)@xmax) onTheGrid = FALSE
																			if (pt2_rotated[1] < extent(hullRaster)@xmin) onTheGrid = FALSE
																			if (pt2_rotated[2] > extent(hullRaster)@ymax) onTheGrid = FALSE
																			if (pt2_rotated[2] < extent(hullRaster)@ymin) onTheGrid = FALSE
																			if (onTheGrid == TRUE)
																				{
																					NAarea = FALSE
																					if (is.na(raster::extract(hullRaster,cbind(pt2_rotated[1],pt2_rotated[2]))))
																						{
																							NAarea = TRUE; twoPointsOnTheGrid = TRUE
																						}
																					if (NAarea == FALSE)
																						{
																							pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																						}
																				}
																			if (counter >= 100)
																				{
																					pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																					pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
																					pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
																				}
																		}
																}
															if (counter < 100)
																{
																	pt2_rotated_list[[g]] = pt2_rotated
																}	else	{
																	pt2_rotated_list[[g]] = pt2
																}
															envValues1[[g]] = raster::extract(envVariable, cbind(pt2_rotated[1],pt2_rotated[2]))
														}
													if (is.null(repulsion))
														{
															pt2 = pt2_rotated_list[[1]]
														}	else	{
															showingOriginalTree = FALSE; showingRotations = FALSE
															if (showingOriginalTree)
																{
																	for (g in 1:dim(fromCoor[[t]])[1])
																		{
																			segments(fromCoor[[t]][g,1], fromCoor[[t]][g,2], toCoor[[t]][g,1], toCoor[[t]][g,2], col="gray30", lwd=0.4)
																		}
																}
															if (showingRotations)
																{
																	for (g in 1:length(pt2_rotated_list))
																		{
																			segments(pt1[1], pt1[2], pt2_rotated_list[[g]][1], pt2_rotated_list[[g]][2], col="gray30", lwd=0.2)
																		}
																}
															envValues2 = unlist(envValues1)
															
															
															envValues2 = unlist(envValues1)
															if (repulsion == TRUE)
																{
																	buffer = envValues2-min(envValues2,na.rm=T)
																	envValues2 = max(envValues2,na.rm=T)-buffer
																}
															vS = envValues2-min(envValues2,na.rm=T)
															probas = vS/sum(vS,na.rm=T)
															probas[which(is.na(probas))] = 0
															if (sum(probas) != 0)
																{
																	index = sample(1:length(envValues2), 1, prob=probas)
																	pt2 = pt2_rotated_list[[index]]
																}
														}
													if (showingPlots == TRUE)
														{
															segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.4)
															points(pt2[1], pt2[2], pch=16, col="black", cex=0.5)
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
											write.csv(temp, paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
											if (showingPlots == TRUE) dev.copy2pdf(file=paste0(randomisationDirectory,"/TreeRandomisation_",t,".pdf"))
										}
									t
								}
						}			
					if (branchRandomisation1 == TRUE)
						{
							fromCoorRand = fromCoor; toCoorRand = toCoor; buffer = list()
							# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
							for (t in 1:nberOfExtractionFiles) {
									newRandomisation = TRUE
									if (file.exists(paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv")))
										{
											newRandomisation = FALSE
											if (overwrite) newRandomisation = TRUE
										}
									if (newRandomisation)
										{
											cat("Randomising tree ",t,"\n",sep="")
											if (showingPlots == TRUE)
												{
													if (t == 1) plotRaster(envVariable, addLegend=T, new=T)
													if (t >= 2) plotRaster(envVariable, addLegend=T, new=T)
													lines(points[hull,], lwd=1.5, col=rgb(222,67,39,220,maxColorValue=255)) # red
													text1 = paste("Randomisation of branch positions, sampled tree ",t,sep="")
													mtext(text1, col="gray30", cex=0.9, line=-1)
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
																	xTranslation = runif(1)*(extent(hullRaster)@xmax-extent(hullRaster)@xmin)
																	yTranslation = runif(1)*(extent(hullRaster)@ymax-extent(hullRaster)@ymin)
																	pt1_translated[1] = pt1[1]+xTranslation; pt1_translated[2] = pt1[2]+yTranslation
																	if (pt1_translated[1] > extent(hullRaster)@xmax)
																		{
																			pt1_translated[1] = pt1_translated[1]-(extent(hullRaster)@xmax-extent(hullRaster)@xmin)
																			xTranslation = xTranslation-(extent(hullRaster)@xmax-extent(hullRaster)@xmin)
																		}
																	if (pt1_translated[2] > extent(hullRaster)@ymax)
																		{
																			pt1_translated[2] = pt1_translated[2]-(extent(hullRaster)@ymax-extent(hullRaster)@ymin)
																			yTranslation = yTranslation-(extent(hullRaster)@ymax-extent(hullRaster)@ymin)
																		}
																	NAarea = FALSE	
																	if (is.na(raster::extract(hullRaster,cbind(pt1_translated[1],pt1_translated[2]))))
																		{
																			NAarea = TRUE
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
															pt2NAarea = TRUE; counter = 0
															while (pt2NAarea == TRUE)
																{
																	counter = counter+1; onTheGrid = TRUE
																	angle = (2*pi)*runif(1); pt2_rotated = rotation(pt1, pt2, angle)
																	if (pt2_rotated[1] > extent(hullRaster)@xmax) onTheGrid = FALSE
																	if (pt2_rotated[1] < extent(hullRaster)@xmin) onTheGrid = FALSE
																	if (pt2_rotated[2] > extent(hullRaster)@ymax) onTheGrid = FALSE
																	if (pt2_rotated[2] < extent(hullRaster)@ymin) onTheGrid = FALSE
																	if (onTheGrid == TRUE)
																		{
																			NAarea = FALSE
																			if (is.na(raster::extract(hullRaster,cbind(pt2_rotated[1],pt2_rotated[2]))))
																				{
																					NAarea = TRUE; twoPointsOnTheGrid = TRUE
																				}
																			if (NAarea == FALSE)
																				{
																					pol_x = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,1]
																					pol_y = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,2]
																					if (point.in.polygon(pt2_rotated[1], pt2_rotated[2], pol_x, pol_y) == 1)
																						{
																							pt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
																						}
																				}
																		}
																	if (counter > 100)
																		{
																			pt2NAarea = FALSE; twoPointsOnTheGrid = FALSE
																			pt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
																			pt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
																		}
																}
														}
													pt2 = pt2_rotated
													if (showingPlots == TRUE)
														{
															segments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.4)
															points(pt2[1], pt2[2], pch=16, col="black", cex=0.5)
														}
													fromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
													toCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
												}
											temp = datas[[t]]
											temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
											write.csv(temp, paste0(randomisationDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
											cdev.copy2pdf(file=paste0(randomisationDirectory,"/TreeRandomisation_",t,".pdf"))
										}
									t
								}
						}	
				}
		}		
}
