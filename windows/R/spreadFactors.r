spreadFactors = function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), pathModel=1, resistances=list(), avgResistances=list(), fourCells=FALSE, 
			 nberOfRandomisations=0, randomProcedure=3, outputName="", showingPlots=FALSE,OS="Unix", juliaCSImplementation=FALSE, 
			 simulations=FALSE, randomisations=FALSE, hull_polygons=list(), onlyTipBranches=FALSE, GLM=FALSE, alternativeQstat=FALSE) {

# variogramModels = list(); hull_polygons = list(); simulations = FALSE; randomisations=FALSE; onlyTipBranches = FALSE; GLM = FALSE; alternativeQstat = FALSE
CA = FALSE
onlyTipBranches = FALSE
impactOnVelocity = TRUE
impactOnDirection = FALSE
if ((pathModel == -1)|(pathModel == 0))
	{
		impactOnVelocity = FALSE
		impactOnDirection = TRUE
	}
# registerDoMC(cores=nberOfCores)
nberOfCores_CS = 1
plottingHistograms = FALSE
commonalityAnalysis = FALSE
if (CA == TRUE) commonalityAnalysis = TRUE
all = FALSE # boolean variable specifying if LR coefficients have be tested with randomisation on all rasters (TRUE) or with one environmental 
	# raster at a time (FALSE). In the case of torus translations, if TRUE, all the rasters will be uniformly translated at each iteration
thetaValue = 1 # numeric value [0 < theta < 20], the degree from which the path randomly deviates from the shortest path
envVariableToRandomise = -1 # index of the environmental data to randomise. By default (= -1), all the data will be randomised
dispersalTimeBoolean = TRUE
date1 = base::date()
preLogTransformation = function(x, m)
	{
		x = 1+(x-m)
	}
logTransformation = function(x, m)
	{
		x = log10(x)
	}
zTransformation = function(x)
	{ 
		x = (x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T))
	}
featureScaling = function(x)
	{ 
		x = (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
	}	
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
newResistances = c(TRUE)
newAvgResistances = c(TRUE)
for (h in 1:length(envVariables0))
	{
		newEnvVariables[[h+1]] = envVariables0[[h]]
		newEnvVariables[[h+1]][newEnvVariables[[h+1]]<0] = NA
		if(resistances[[h]] == TRUE)
			{
				fric = "_R"
			}	else		{
				fric = "_C"
			}
		names(newEnvVariables[[h+1]]) = paste(names(newEnvVariables[[h+1]]), fric, sep="")
		if (length(resistances) > 0)
			{
				newResistances = c(newResistances, resistances[[h]])
				newAvgResistances = c(newAvgResistances, avgResistances[[h]])
			}
	}
envVariables = newEnvVariables
resistances = newResistances
avgResistances = newAvgResistances
distanceMatrix = FALSE
straightLineDistance = FALSE
leastCostDistance = FALSE
rSPDistance = FALSE
commuteDistance = FALSE
randomWalkDistance = FALSE
torusRandomisations = FALSE
rastersSimulations = FALSE
externalRandomisations = FALSE
externalSimulations = FALSE
branchRandomisation3 = FALSE
branchRandomisation2 = FALSE
branchRandomisation1 = FALSE
distPermutations = FALSE
if (pathModel == 1) straightLineDistance = TRUE
if (pathModel == 2) leastCostDistance = TRUE
if (pathModel == 3) randomWalkDistance = TRUE
if (pathModel == 4) commuteDistance = TRUE
if (pathModel == 5) rSPDistance = TRUE
if (randomProcedure == 1) externalRandomisations = TRUE
if (randomProcedure == 2) externalSimulations = TRUE
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
if (randomProcedure == 7) distPermutations = TRUE
nberOfConnections = matrix(nrow=1, ncol=nberOfExtractionFiles)	
totalnberOfConnections = 0
if ((simulations == FALSE)&(randomisations == FALSE))
	{
		extractionFileName = "TreeExtractions"
	}
if ((simulations == TRUE)&(randomisations == FALSE))
	{
		extractionFileName = "TreeSimulations"
	}
if ((simulations == FALSE)&(randomisations == TRUE))
	{
		extractionFileName = "TreeRandomisation"
	}
if (nchar(localTreesDirectory) == 0)
	{
		data = read.csv(paste(extractionFileName,"_1.csv",sep=""), header=T, dec=".")
	}	else		{
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
if (impactOnVelocity == TRUE)
	{
		distances = list()
	}
if (impactOnDirection == TRUE)
	{
		meanEnvValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		rateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		meanDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
	}
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
		if (onlyTipBranches == TRUE)
			{
				indices = which(!data[,"node2"]%in%data[,"node1"])
				data = data[indices,]
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
		if (impactOnVelocity == TRUE)
			{
				distances[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
			}
		if (("treeID"%in%colnames(data)) == TRUE)
			{
				treeIDs[[t]] = data[1,"treeID"]
			}	else	{
				treeIDs[[t]] = "noTreeID"
			}	
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
		hullRasters[[h]] = raster::mask(hullRasters[[h]], sps, snap="out")
		hullRasters[[h]][!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
		names(hullRasters[[h]]) = gsub(".asc","",names(envVariables[[h]]))
		names(hullRasters[[h]]) = gsub(".tif","",names(envVariables[[h]]))
		names(hullRasters[[h]]) = gsub(".gri","",names(envVariables[[h]]))
	}
if (randomWalkDistance == TRUE)
	{
		extensions = rep("", length(envVariables))
		if ("CS_rasters"%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), "CS_rasters"))
		for (h in 1:length(envVariables))
			{
				extensions[h] = ".asc"
				if (round(res(hullRasters[[h]])[1],10) != round(res(hullRasters[[h]])[2],10)) extensions[h] = ".tif"
				name = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
				writeRaster(hullRasters[[h]], name, overwrite=T)
			}
	}

## 1. Estimation of the different environmental distances or differences associated with each branch

if (distanceMatrix == F)
	{
		for (h in 1:length(envVariables))
			{
				if (impactOnVelocity == TRUE)
					{
						if (straightLineDistance == TRUE) cat("Computing environmental distances (straight-line path model) for ", names(envVariables[[h]])[1], "\n", sep="")
						if (leastCostDistance == TRUE) cat("Computing environmental distances (least-cost path model) for ", names(envVariables[[h]])[1], "\n", sep="")
						if (randomWalkDistance == TRUE) cat("Computing environmental distances (Circuitscape path model) for ", names(envVariables[[h]])[1], "\n", sep="")	
						if (fourCells == TRUE) directions = 4
						if (fourCells == FALSE) directions = 8
						if (resistances[h] == FALSE)
							{
								if (leastCostDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], mean, directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
									}
								if (commuteDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], mean, directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
									}
								if (rSPDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], mean, directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
									}
							}	else	{
								if (leastCostDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
									}
								if (rSPDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
									}
								if (randomWalkDistance == TRUE)
									{
										trEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
										trEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
									}	
							}
						if (showingPlots == FALSE)
							{
								buffer = list()
								# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
								for (t in 1:nberOfExtractionFiles) {
										mat = matrix(nrow=nberOfConnections[t], ncol=1)
										if (straightLineDistance == TRUE)
											{
												linesList = list()
												for (i in 1:length(fromCoor[[t]][,1]))
													{
														points = rbind(fromCoor[[t]][i,], toCoor[[t]][i,])
														linesList[[i]] = Lines(list(Line(points)),i)
													}
												lines = SpatialLines(linesList)
												extractions = raster::extract(hullRasters[[h]], lines)
												for (i in 1:length(fromCoor[[t]][,1]))
													{
														mat[i] = sum(extractions[[i]], na.rm=T)
													}
											}
										if (leastCostDistance == TRUE)
											{
												mat[] = diag(costDistance(trEnvVariableCorr, fromCoor[[t]], toCoor[[t]]))									
											}
										if (randomWalkDistance == TRUE)
											{
												envVariableName = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
												branchesNotNA = which(!((is.na(raster::extract(hullRasters[[h]], fromCoor[[t]][])))
																	  |(is.na(raster::extract(hullRasters[[h]], toCoor[[t]][])))))
												fromCoor_temp = fromCoor[[t]][branchesNotNA,]; toCoor_temp = toCoor[[t]][branchesNotNA,]
												if (juliaCSImplementation == FALSE)
													{
														mat[branchesNotNA,1] = circuitScape1(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
																							 fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
														# if (-777%in%mat[])
															# {
																# mat[branchesNotNA,1] = circuitScape1(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
																									 # fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
															# }
													}	else	{
														mat[branchesNotNA,1] = circuitScape2(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
																							 fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
													}
											}	
										if (commuteDistance == TRUE)
											{
												for (i in 1:length(fromCoor[[t]][,1]))
													{
														spatialPoints = SpatialPoints(cbind(c(fromCoor[[t]][i,1], toCoor[[t]][i,1]), c(fromCoor[[t]][i,2], toCoor[[t]][i,2])))
														mat[i] = commuteDistance(trEnvVariableCorr, spatialPoints)
													}	
											}
										if (rSPDistance == TRUE)
											{
												mat[] = diag(rSPDistance(trEnvVariableCorr, fromCoor[[t]], toCoor[[t]], theta=thetaValue, totalNet="total", method=1))
											}		
										colnames(mat) = names(envVariables[[h]])
										buffer[[t]] = mat
										# mat
									}		
								for (t in 1:length(buffer))
									{
										buffer[[t]][!is.finite(buffer[[t]][])] = NA # least-cost case
										buffer[[t]][buffer[[t]][]==-1] = NA # CircuitScape (RW) case
										buffer[[t]][buffer[[t]][]==-777] = NA # NA value in Circuitscape
										distances[[t]][,h] = buffer[[t]][]
									}		
							}	else	{
								for (t in 1:nberOfExtractionFiles)
									{
										plotRaster(envVariables[[h]], addLegend=T)
										lines(points[hull,], lwd=0.5, col="black")
										text1 = paste(names(envVariables[[h]])[1],", sampled tree ",t,sep="")
										if (straightLineDistance == TRUE) text2 = "Computing straight-line distances..."
										if (leastCostDistance == TRUE) text2 = "Computing least-cost distances..."
										if ((randomWalkDistance == TRUE)|(commuteDistance == TRUE)) text2 = "Computing Circuitscape distances..."
										if (rSPDistance == TRUE) text2 = "Computing randomised shortest path distances..."
										mtext(text1, col="black", cex=0.7, line=0)
										mtext(text2, col="red", cex=0.7, line=-1)
										for (j in 1:nberOfConnections[t])
											{
												if (j == 1)
													{
														points(fromCoor[[t]][j,1], fromCoor[[t]][j,2], pch=16, col="black", cex=0.5)
													}
												segments(fromCoor[[t]][j,1], fromCoor[[t]][j,2], toCoor[[t]][j,1], toCoor[[t]][j,2], col="black", lwd=0.3)
												points(toCoor[[t]][j,1], toCoor[[t]][j,2], pch=16, col="black", cex=0.5)
												if (straightLineDistance == TRUE)
													{	
														linesList = list()
														points = rbind(fromCoor[[t]][j,], toCoor[[t]][j,])
														linesList[[1]] = Lines(list(Line(points)),j)
														lines = SpatialLines(linesList)
														extractions = raster::extract(envVariables[[h]], lines)
														distances[[t]][j,h] = sum(extractions[[1]], na.rm=T)		
													}
												if (leastCostDistance == TRUE)
													{
														fromCoorJ = cbind(fromCoor[[t]][j,1], fromCoor[[t]][j,2])
														toCoorJ = cbind(toCoor[[t]][j,1], toCoor[[t]][j,2])
														distances[[t]][j,h] = costDistance(trEnvVariableCorr, fromCoorJ, toCoorJ)
													}
												if (randomWalkDistance == TRUE)
													{
														envVariableName = paste("CS_rasters/",names(envVariables[[h]]),"_",outputName,"_cs",extensions[h],sep="")
														fromC = matrix(nrow=1, ncol=2); toC = matrix(nrow=1, ncol=2)
														fromC[] = fromCoor[[t]][j,]; toC[] = toCoor[[t]][j,]
														if (juliaCSImplementation == FALSE)
															{
																distances[[t]][j,h] = circuitScape1(envVariables[[h]], envVariableName, resistances[[h]], avgResistances[[h]], 
																						   			fourCells, fromC, toC, OS, outputName,t, nberOfCores_CS)
															}	else		{
																distances[[t]][j,h] = circuitScape2(envVariables[[h]], envVariableName, resistances[[h]], avgResistances[[h]], 
																						   			fourCells, fromC, toC, OS, outputName,t, nberOfCores_CS)
															}
													}
												if (commuteDistance == TRUE)
													{
														spatialPoints = SpatialPoints(cbind(c(fromCoor[[t]][j,1], toCoor[[t]][j,1]), c(fromCoor[[t]][j,2], toCoor[[t]][j,2])))
														distances[[t]][j,h] = commuteDistance(trEnvVariableCorr, spatialPoints)
													}
												if (rSPDistance == TRUE)
													{
														distances[[t]][j,h] = rSPDistance(trEnvVariableCorr, fromCoor[[t]][j,], toCoor[[t]][j,],
																						  theta=thetaValue, totalNet="total", method=1)
													}	
												if (is.na(as.numeric(distances[[t]][j,h])) == FALSE)
													{
														if (!is.finite(distances[[t]][j,h]))
															{
																distances[[t]][j,h] = NA
															}	else	{
																if (distances[[t]][j,h] == -1) distances[[t]][j,h] = NA	
															}
													}
											}
										dev.off()
									}
							}
					}
				if (impactOnDirection == TRUE)
					{
						for (t in 1:nberOfExtractionFiles)
							{
								envValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
								for (i in 1:length(ancestralNodes))
									{
										ancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
										envValues = envValues + raster::extract(envVariables[[h]], cbind(fromCoor[[t]][ancestralBranch,1],fromCoor[[t]][ancestralBranch,2]))
									}
								meanEnvValues[t,h] = envValues + mean(raster::extract(envVariables[[h]], toCoor[[t]]), na.rm=T)
								diffs = raster::extract(envVariables[[h]], fromCoor[[t]])-raster::extract(envVariables[[h]], toCoor[[t]])
								rateOfPositiveDifferences[t,h] = sum(diffs[!is.na(diffs)] > 0)/length(diffs[!is.na(diffs)])
								meanDifferences[t,h] = mean(diffs, na.rm=T)
							}
					}
			}					
	}

envVariableNames = names(envVariables[[1]])	
for (h in 2:length(envVariables))
	{
		envVariableNames = cbind(envVariableNames, names(envVariables[[h]]))
	}
for (t in 1:nberOfExtractionFiles)
	{
		if (impactOnVelocity == TRUE) colnames(distances[[t]]) = envVariableNames
	}

# To generate scatterplots:
if ((showingPlots == TRUE) & (impactOnVelocity == TRUE))
	{
		# for (h in 1:length(envVariables))
			# {
				# fileName = paste(outputName,"_scatterplot_",names(envVariables[[h]]),".pdf",sep="")
				# dists = distances[[1]][,h]
				# times = dispersalTime[[1]]
				# if (nberOfExtractionFiles > 1)
					# {
						# for (t in 2:nberOfExtractionFiles)
							# {
								# dists = c(dists, distances[[t]][,h])
								# times = c(times, dispersalTime[[t]])
							# }
					# }
				# dev.new(width=5, height=5); par(mar=c(3.5,3.5,4,2), mgp=c(1.6,0.5,0))	
				# formula = paste("times ~ dists", sep="")
				# form = as.formula(formula); LM = lm(form)
				# title = paste("Linear regression R2 = ",round(summary(LM)$r.squared,3),sep="")
				# plot(dists, times, xlab="environmental distances", ylab="dispersal times", main=title, cex.main=0.7, cex.axis=0.6, lwd=0.8, pch=19, cex=0.3, cex.lab=0.7)
				# abline(LM, lwd=1.5, col="red")
				# dev.copy2pdf(file=fileName)
				# dev.off()
			# }
	}
if ((nberOfExtractionFiles == 1) & (impactOnVelocity == TRUE))
	{
		fileName = paste(outputName, "_env_distances.txt", sep="")
		mat = cbind(dispersalTime[[1]][], distances[[1]][,1:length(envVariables)])
		columnNames = "dispersal_times"
		for (h in 1:length(envVariables))
			{
				columnNames = cbind(columnNames, names(envVariables[[h]]))
			}
		colnames(mat) = columnNames
		write.table(mat, file=fileName, row.names=F, quote=F, sep="\t")
	}
if ((file.exists(outputName)) & (impactOnVelocity == TRUE))
	{	
		for (t in 1:nberOfExtractionFiles)
			{
				fileName = paste(outputName, "/", outputName, "_tree", t, "_env_distances.txt", sep="")
				mat = cbind(dispersalTime[[t]][], distances[[t]][,1:length(envVariables)])
				columnNames = cbind("dispersal_times")
				for (h in 1:length(envVariables))
					{
						columnNames = cbind(columnNames, names(envVariables[[h]]))
					}
				colnames(mat) = columnNames
				write.table(mat, file=fileName, row.names=F, quote=F, sep="\t")
			}
	}

## 2. Linear regressions of dispersal time vs. environmental distances

if (impactOnVelocity == TRUE)
	{
		# realUniLRresiduals = list()
		realUniLRcoefficients1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRcoefficients2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquarePValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		if (GLM == TRUE)
			{
				realMultiGLMcoefficients = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				realMultiGLMresiduals = matrix(nrow=length(dispersalTime[[1]]), ncol=nberOfExtractionFiles)
				realMultiGLMcoefficientPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				realMultiGLMCAuniqueContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				realMultiGLMCAcommonContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
			}
		colNames = c(); for (t in 1:nberOfExtractionFiles) colNames = c(colNames, paste("tree_",treeIDs[[t]],sep=""))
		if (GLM == TRUE) colnames(realMultiGLMresiduals) = colNames
		for (h in 1:length(envVariables))
			{
				nRowsMax = length(dispersalTime[[1]])
				if (nberOfExtractionFiles > 1)
					{
						for (t in 2:nberOfExtractionFiles)
							{
								if (nRowsMax < length(dispersalTime[[t]])) nRowsMax = length(dispersalTime[[t]])
							}
					}
				# residuals = matrix(nrow=nRowsMax, ncol=nberOfExtractionFiles)
				for (t in 1:nberOfExtractionFiles)
					{
						distVariables = paste("dispersalTime[[",t,"]]"," ~ distances[[",t,"]][,",h,"]",sep="")
						form = as.formula(distVariables)
						LM = lm(form)
						realUniLRcoefficients1[t,h] = summary(LM)$coefficients[2,"Estimate"]
						realUniLRRsquares1[t,h] = summary(LM)$r.squared
						f = summary(LM)$fstatistic
						if (is.numeric(f))
							{
								p = pf(f[1],f[2],f[3],lower.tail=F)
								attributes(p) = NULL; realUniLRRsquarePValues[t,h] = p
							}	else		{
								realUniLRRsquarePValues[t,h] = NA 
							}
						nonNaNdispersalTimes = dispersalTime[[t]][!is.na(distances[[t]][,1])]
						nonNaNdistances1 = distances[[t]][,1][!is.na(distances[[t]][,1])]
						nonNaNdistances2 = distances[[t]][,h][!is.na(distances[[t]][,1])]
						distVariables = paste("nonNaNdistances2 ~ nonNaNdistances1",sep="")
						LM1 = lm(as.formula(distVariables)); residuals_LM1 = LM1$residuals; indices = which(!is.na(residuals_LM1))
						distVariables = paste("nonNaNdispersalTimes[indices] ~ nonNaNdistances2[indices] + residuals_LM1",sep="")
						LM2 = lm(as.formula(distVariables)); realUniLRRsquares2[t,h] = summary(LM2)$r.squared
						realUniLRcoefficients2[t,h] = summary(LM2)$coefficients[2,"Estimate"]
						# if (length(residuals[,t]) == length(studres(LM)))
							# {
								# residuals[,t] = studres(LM)
							# }	else		{
								# for (i in 1:length(studres(LM)))
									# {
										# residuals[i,t] = studres(LM)[i]
									# }
							# }
					}
				# realUniLRresiduals[[h]] = residuals
			}
		if (GLM == TRUE)
			{
				for (t in 1:nberOfExtractionFiles)
					{
						distVariables = ""
						matMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
						matMultiGLM[,1] = dispersalTime[[t]]
						matMultiGLMNames = c("dispersalTime")
						for (h in 1:length(envVariables))
							{
								matMultiGLM[,h+1] = distances[[t]][,h]
								matMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
								if (h == 1)
									{
										distVariables = paste("dispersalTime ~ ", names(envVariables[[h]]), sep="")
									}	else	{
										distVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
									}
							}
						matMultiGLM = as.data.frame(matMultiGLM)
						names(matMultiGLM) = matMultiGLMNames
						m = min(unlist(lapply(as.matrix(matMultiGLM),min)), na.rm=T)	
						matMultiGLM = lapply(matMultiGLM, preLogTransformation, m)
						matMultiGLM = lapply(matMultiGLM, logTransformation)
						matMultiGLM = lapply(matMultiGLM, zTransformation)
						matMultiGLM = as.data.frame(matMultiGLM)
						names(matMultiGLM) = matMultiGLMNames
						form = as.formula(distVariables)
						multiGLM = stats::glm(form, data=matMultiGLM)
						if (length(realMultiGLMresiduals[,t]) == length(studres(multiGLM)))
							{
								realMultiGLMresiduals[,t] = studres(multiGLM)
							}	else	{
								for (i in 1:length(studres(multiGLM)))
									{
										realMultiGLMresiduals[i,t] = studres(multiGLM)[i]
									}
							}
						if (commonalityAnalysis == TRUE)
							{
								CA = calc.yhat(multiGLM, prec=5)	
							}
						for (h in 1:length(envVariables))
							{
								names(multiGLM$coefficients)[1+h] = names(envVariables[[h]])
							}
						for (h in 1:length(envVariables))
							{	
								name = names(envVariables[[h]]); nameInside = FALSE
								for (i in 1:length(summary(multiGLM)[]$coefficients[,"Estimate"]))
									{
										if (name == rownames(summary(multiGLM)[]$coefficients)[i])
											{
												 nameInside = TRUE
											}
									}
								if (nameInside == TRUE)
									{	
										realMultiGLMcoefficients[t,h] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
										realMultiGLMcoefficientPValues[t,h] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
										if (commonalityAnalysis == TRUE)
											{
												realMultiGLMCAuniqueContributions[t,h] = CA$PredictorMetrics[h,"Unique"]
												realMultiGLMCAcommonContributions[t,h] = CA$PredictorMetrics[h,"Common"]	
											}
									}	
							}
					}
			}
				
		uniLRmedianRsquares = matrix(nrow=1, ncol=length(envVariables))
		uniLRmedianRsquarePValues = matrix(nrow=1, ncol=length(envVariables))
		uniLRmedianDeltaRsquares = matrix(nrow=1, ncol=length(envVariables))
		if (GLM == TRUE)
			{
				multiGLMmedianCoefficients = matrix(nrow=1, ncol=length(envVariables))
				multiGLMmedianCoefficientPValues = matrix(nrow=1, ncol=length(envVariables))
				multiGLMmedianCAuniqueContributions = matrix(nrow=1, ncol=length(envVariables))
				multiGLMmedianCAcommonContributions = matrix(nrow=1, ncol=length(envVariables))
			}
		for (h in 1:length(envVariables))
			{
				uniLRmedianRsquares[1,h] = median(realUniLRRsquares1[,h], na.rm=T)
				uniLRmedianRsquarePValues[1,h] = median(realUniLRRsquarePValues[,h], na.rm=T)
				uniLRmedianDeltaRsquares[1,h] = median(realUniLRRsquares1[,h]-realUniLRRsquares1[,1], na.rm=T)
				if (GLM == TRUE)
					{
						multiGLMmedianCoefficients[1,h] = median(realMultiGLMcoefficients[,h], na.rm=T)
						multiGLMmedianCoefficientPValues[1,h] = median(realMultiGLMcoefficientPValues[,h], na.rm=T)
						if (commonalityAnalysis == TRUE)
							{
								multiGLMmedianCAuniqueContributions[1,h] = median(realMultiGLMCAuniqueContributions[,h], na.rm=T)
								multiGLMmedianCAcommonContributions[1,h] = median(realMultiGLMCAcommonContributions[,h], na.rm=T)
							}
					}
				if (h != 1)
					{
						realUniDeltaRsquares1[,h] = (realUniLRRsquares1[,h]-realUniLRRsquares1[,1])
						realUniDeltaRsquares2[,h] = (realUniLRRsquares2[,h]-realUniLRRsquares1[,1])
					}
			}
			
		if ((plottingHistograms == TRUE)&(nberOfExtractionFiles > 1))
			{
				if (GLM == TRUE)
					{
						if (commonalityAnalysis == TRUE)
							{
								fileName1 = paste(outputName, "_LR-GLM-CA_results.pdf", sep="")
								pdf(fileName1, width=(4*(length(envVariables))), height=(5*4)); par(mfrow=c(6,(length(envVariables))))
							}	else	{
								fileName1 = paste(outputName, "_LR-GLM_results.pdf", sep="")
								pdf(fileName1, width=(4*(length(envVariables))), height=(3*4)); par(mfrow=c(4,(length(envVariables))))
							}
					}	else	{
						fileName1 = paste(outputName, "_linear_regression_results.pdf", sep="")
						pdf(fileName1, width=(4*(length(envVariables))), height=(2.25*4)); par(mfrow=c(3,(length(envVariables))))
					}
				breakList_1 = (0:50)/50; breakList_2 = (0:50)/100
				for (h in 1:length(envVariables))
					{
						xMin = min(realUniLRcoefficients1[,h]); xMax = max(realUniLRcoefficients1[,h]) # breaks=seq(min,max,by=(max-min)/50)
						hist(realUniLRcoefficients1[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50), 
							 xlab="Univariate LR coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
					}
				xMin = min(realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,1], na.rm=T)	
				for (h in 2:length(envVariables))
					{
						if (xMin > min(realUniLRRsquares1[,h], na.rm=T))
							{
								xMin = min(realUniLRRsquares1[,h], na.rm=T)
							}
						if (xMax < max(realUniLRRsquares1[,h], na.rm=T))
							{
								xMax = max(realUniLRRsquares1[,h], na.rm=T)
							}
					}
				for (h in 1:length(envVariables))
					{
						min = min(realUniLRRsquares1[,h]); max = max(realUniLRRsquares1[,h]) # breaks=seq(min,max,by=(max-min)/50)
						hist(realUniLRRsquares1[,h], freq=T, xlim=c(0,xMax), breaks=seq(0,xMax,by=(xMax-0)/50), xlab="Univariate LR R2's", 
							 main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
					}
				plot.new()
				xMin = min(realUniLRRsquares1[,2]-realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,2]-realUniLRRsquares1[,1], na.rm=T)
				for (h in 2:length(envVariables))
					{
						if (xMin > min(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T))
							{
								xMin = min(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T)
							}
						if (xMax < max(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T))
							{
								xMax = max(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T)
							}
					}
				for (h in 2:length(envVariables))
					{
						hist(realUniDeltaRsquares1[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50), xlab="Univariate LR delta R2 (Q)", 
							 main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
					}
				if (GLM == TRUE)
					{
						xMin = min(realMultiGLMcoefficients[,1], na.rm=T); xMax = max(realMultiGLMcoefficients[,1], na.rm=T)
						for (h in 2:length(envVariables))
							{
								if (xMin > min(realMultiGLMcoefficients[,h], na.rm=T))
									{
										xMin = min(realMultiGLMcoefficients[,h], na.rm=T)
									}
								if (xMax < max(realMultiGLMcoefficients[,h], na.rm=T))
									{
										xMax = max(realMultiGLMcoefficients[,h], na.rm=T)
									}
							}
						for (h in 1:length(envVariables))
							{
								if (is.na(mean(realMultiGLMcoefficients[,h])))
									{	
										plot.new()
									}	else		{
										hist(realMultiGLMcoefficients[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
											 xlab="Multivariate GLM coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
									}
							}	
						if (commonalityAnalysis == TRUE)
							{
								xMin = min(c(realMultiGLMCAuniqueContributions[,1],realMultiGLMCAcommonContributions[,1]), na.rm=T)
								xMax = max(c(realMultiGLMCAuniqueContributions[,1],realMultiGLMCAcommonContributions[,1]), na.rm=T)
								for (h in 2:length(envVariables))
									{
										if (xMin > min(realMultiGLMCAuniqueContributions[,h], na.rm=T))
											{
												xMin = min(realMultiGLMCAuniqueContributions[,h], na.rm=T)
											}
										if (xMin > min(realMultiGLMCAcommonContributions[,h], na.rm=T))
											{
												xMin = min(realMultiGLMCAcommonContributions[,h], na.rm=T)
											}	
										if (xMax < max(realMultiGLMCAuniqueContributions[,h], na.rm=T))
											{
												xMax = max(realMultiGLMCAuniqueContributions[,h], na.rm=T)
											}
										if (xMax < max(realMultiGLMCAcommonContributions[,h], na.rm=T))
											{
												xMax = max(realMultiGLMCAcommonContributions[,h], na.rm=T)
											}	
									}	
								for (h in 1:length(envVariables))
									{
										if (is.na(mean(realMultiGLMcoefficients[,h])))
											{	
												plot.new()
											}	else		{
												hist(realMultiGLMCAuniqueContributions[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
												xlab="Multivariate GLM-CA unique contributions", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}	
									}
								for (h in 1:length(envVariables))
									{
										if (is.na(mean(realMultiGLMcoefficients[,h])))
											{	
												plot.new()
											}	else	{
												hist(realMultiGLMCAcommonContributions[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
												xlab="Multivariate GLM-CA common contributions", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}	
									}
							}
					}
				dev.off()
			}
		if (nberOfExtractionFiles > 1)
			{
				if (GLM == TRUE)
					{
						if (commonalityAnalysis == TRUE)
							{
								fileName2 = paste(outputName,"_LR-GLM-CA_results.txt",sep="")
							}	else	{
								fileName2 = paste(outputName,"_LR-GLM_results.txt",sep="")
							}
					}	else	{
						fileName2 = paste(outputName, "_linear_regression_results.txt", sep="")
					}
				uniLRcoefficientsNames1 = c(); uniLRRsquaresNames1 = c(); uniLRdeltaRsquaresNames1 = c(); uniLRpValuesNames = c()
				uniLRcoefficientsNames2 = c(); uniLRRsquaresNames2 = c(); uniLRdeltaRsquaresNames2 = c()
				if (GLM == TRUE)
					{
						multiGLMcoefficientsNames = c(); multiGLMpValuesNames = c()
						mat = cbind(realUniLRcoefficients1[,1:dim(realUniLRcoefficients1)[2]],realUniLRRsquares1[,1:dim(realUniLRRsquares1)[2]],
									realUniDeltaRsquares1[,2:dim(realUniDeltaRsquares1)[2]],realMultiGLMcoefficients[,1:dim(realMultiGLMcoefficients)[2]])
						if (commonalityAnalysis == TRUE)
							{
								multiGLMCAuniqueContributionsNames = c(); multiGLMCAcommonContributionsNames = c()
								mat = cbind(mat, realMultiGLMCAuniqueContributions[,1:dim(realMultiGLMCAuniqueContributions)[2]],
												 realMultiGLMCAcommonContributions[,1:dim(realMultiGLMCAcommonContributions)[2]])
							}
					}	else	{
						mat = cbind(realUniLRcoefficients1[,1:dim(realUniLRcoefficients1)[2]], realUniLRRsquares1[,1:dim(realUniLRRsquares1)[2]], 
																	   realUniDeltaRsquares1[,2:dim(realUniDeltaRsquares1)[2]])
						if (alternativeQstat == TRUE)
							{
								mat = cbind(mat, realUniLRcoefficients2[,1:dim(realUniLRcoefficients2)[2]], realUniLRRsquares2[,1:dim(realUniLRRsquares2)[2]], 
																	   realUniDeltaRsquares2[,2:dim(realUniDeltaRsquares2)[2]])
							}
					}
				for (h in 1:length(envVariables))
					{
						if (alternativeQstat == TRUE)
							{
								uniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("Univariate_LR1_coefficients_", names(envVariables[[h]]), sep=""))	
								uniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("Univariate_LR1_R2_", names(envVariables[[h]]), sep=""))				
								uniLRcoefficientsNames2 = cbind(uniLRcoefficientsNames2, paste("Univariate_LR2_coefficients_", names(envVariables[[h]]), sep=""))	
								uniLRRsquaresNames2 = cbind(uniLRRsquaresNames2, paste("Univariate_LR2_R2_", names(envVariables[[h]]), sep=""))	
								if (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1, paste("Univariate_LR1_delta_R2_", names(envVariables[[h]]), sep=""))
								if (h > 1) uniLRdeltaRsquaresNames2 = cbind(uniLRdeltaRsquaresNames2, paste("Univariate_LR2_delta_R2_", names(envVariables[[h]]), sep=""))
							}	else		{
								uniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("Univariate_LR_coefficients_", names(envVariables[[h]]), sep=""))	
								uniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("Univariate_LR_R2_", names(envVariables[[h]]), sep=""))				
								if (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1, paste("Univariate_LR_delta_R2_", names(envVariables[[h]]), sep=""))
							}
						if (GLM == TRUE)
							{
								multiGLMcoefficientsNames = cbind(multiGLMcoefficientsNames, paste("Multivariate_GLM_coefficients_", names(envVariables[[h]]), sep=""))
								if (commonalityAnalysis == TRUE)
									{
										multiGLMCAuniqueContributionsNames = cbind(multiGLMCAuniqueContributionsNames,
											paste0("Multivariate_GLM-CA_unique_contributions_",names(envVariables[[h]])))
										multiGLMCAcommonContributionsNames = cbind(multiGLMCAcommonContributionsNames,
											paste0("Multivariate_GLM-CA_common_contributions_",names(envVariables[[h]])))
									}
							}
					}
				if (GLM == TRUE)
					{
						names = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRpValuesNames, uniLRdeltaRsquaresNames1, multiGLMcoefficientsNames, multiGLMpValuesNames)
						if (commonalityAnalysis == TRUE)
							{
								names = cbind(names, multiGLMCAuniqueContributionsNames, multiGLMCAcommonContributionsNames)
							}
					}	else	{
						if (alternativeQstat == TRUE)
							{
								names = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1, 
											  uniLRcoefficientsNames2, uniLRRsquaresNames2, uniLRdeltaRsquaresNames2)
							}	else		{
								names = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1)
							}
					}
				colnames(mat) = names; write.table(mat, file=fileName2, row.names=F, quote=F, sep="\t")
			}
	}

## 3. randomisation step (permutations of least-cost resistances, torus translation or simulation of environmental rasters)

if (nberOfRandomisations > 0)
	{	
		if (impactOnVelocity == TRUE)
			{
				uniLRRsquaresLower1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquaresHigher1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquaresRandomisationPValues1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresLower1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresHigher1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresRandomisationPValues1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresRandomisationBFs1 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)	
				uniLRRsquarePValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquarePValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquarePValuesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquaresLower2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquaresHigher2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRRsquaresRandomisationPValues2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresLower2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresHigher2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresRandomisationPValues2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				uniLRdeltaRsquaresRandomisationBFs2 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)	
				if (GLM == TRUE)
					{
						multiGLMcoefficientsLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMcoefficientsHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMcoefficientsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMcoefficientPValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMcoefficientPValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMcoefficientPValuesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAuniqueContributionsLower = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAuniqueContributionsHigher = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAuniqueContributionsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAcommonContributionsLower = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAcommonContributionsHigher = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						multiGLMCAcommonContributionsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
					}
				# Matrices for e.g. the "median" p-value tests (number of lines = number of randomisations):		
				uniLRmedianRsquaresSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
				uniLRmedianRsquarePValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
				uniLRmedianDeltaRsquaresSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
				# uniLRdeltaRsquaresPValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
				if (GLM == TRUE)
					{
						multiGLMmedianCoefficientsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
						multiGLMmedianCoefficientPValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
						multiGLMmedianCAuniqueContributionsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
						multiGLMmedianCAcommonContributionsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
					}
				if (GLM == TRUE)
					{
						for (i in 1:dim(realMultiGLMcoefficients)[2])
							{
								for (j in 1:dim(realMultiGLMcoefficients)[1])
									{
										if (is.na(realMultiGLMcoefficients[j,i]))
											{
												multiGLMcoefficientsLower[j,i] = NA; multiGLMcoefficientsHigher[j,i] = NA
												multiGLMcoefficientsRandomisationPValues[j,i] = NA
												multiGLMcoefficientPValuesLower[j,i] = NA; multiGLMcoefficientPValuesHigher[j,i] = NA; 
												multiGLMcoefficientPValuesRandomisationPValues[j,i] = NA
											}
									}
							}
					}
			}
		if (impactOnDirection == TRUE)
			{
				meanEnvValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanEnvValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanEnvValuesRandomisationPValues = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanEnvValuesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
				rateOfPositiveDifferencesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				rateOfPositiveDifferencesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				rateOfPositiveDifferencesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				rateOfPositiveDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
				meanDifferencesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanDifferencesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanDifferencesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
				meanDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)	
			}
		for (s in 1:nberOfRandomisations)
			{
				if (impactOnVelocity == TRUE)
					{
						uniLRdeltaRsquaresRand = matrix(0, nrow=nberOfExtractionFiles, ncol=nberOfRandomisations)
						distancesSim = list()
						for (t in 1:nberOfExtractionFiles)
							{
								distancesSim[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
							}
					}
				# 3.1. randomisation step:
				if ((externalRandomisations == TRUE)|(externalSimulations == TRUE))
					{
						simRasters = list(); simRasters = hullRasters
						if (externalRandomisations == TRUE)
							{
								cat("Analysis on external randomisation ",s,"\n",sep="")
								extractionFileName = "TreeRandomisation"
							}
						if (externalSimulations == TRUE)
							{
								cat("Analysis on external simulation ",s,"\n",sep="")
								extractionFileName = "TreeSimulations"
							}
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
								if (onlyTipBranches == TRUE)
									{
										indices = which(!data[,"node2"]%in%data[,"node1"])
										data = data[indices,]
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
								# if (impactOnVelocity == TRUE) # ERROR ??
									# {
										# distances[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
									# }
								if (("treeID"%in%colnames(data)) == TRUE)
									{
										treeIDs[[t]] = data[1,"treeID"]
									}	else	{
										treeIDs[[t]] = "noTreeID"
									}	
							}
						fromCoorRand = fromCoor; toCoorRand = toCoor	
					}
				if (torusRandomisations == TRUE)
					{
						simRasters = list()
						cat("Analysis on torus translation randomisation ",s,"\n",sep="")	
						for (h in 2:length(hullRasters))
							{
								simRasters[[h]] = torusRandomisation(hullRasters[[h]])
								if (showingPlots == TRUE)
									{
										plotRaster(merge(simRasters[[h]],envVariables[[h]]), addLegend=T)
										lines(points[hull,], lwd=0.5, col="black")
										text1 = paste("torus randomisation ",s," of ",names(hullRasters[[h]])[1],sep="")
										mtext(text1, col="black", cex=0.7, line=0)
									}
							}
						fromCoorRand = fromCoor; toCoorRand = toCoor						
					}
				if (rastersSimulations == TRUE)
					{	
						simRasters = list()
						cat("Analysis on raster simulation ", s, "\n", sep="")
						for (h in 2:length(hullRasters))
							{
								simRasters[[h]] = rasterSimulation(hullRasters[[h]], variogramModels[[h-1]]) 
								if (showingPlots == TRUE)
									{
										plotRaster(merge(simRasters[[h]],envVariables[[h]]), addLegend=T)
										lines(points[hull,], lwd=0.5, col="black")
										text1 = paste("raster simulation ",s," for ",names(hullRasters[[h]])[1],sep="")
										mtext(text1, col="black", cex=0.7, line=0)
									}
							}
						fromCoorRand = fromCoor; toCoorRand = toCoor
					}
				if (branchRandomisation3 == TRUE)
					{
						simRasters = list()
						simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
							{
								counter1 = 0
								twoPointsOnTheGrid = FALSE
								while (twoPointsOnTheGrid == FALSE)
									{
										twoPointsOnTheGrid = TRUE; counter1 = counter1+1
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
																									}	else		{
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
							}
					}
				if (branchRandomisation2 == TRUE)
					{
						simRasters = list()
						simRasters = hullRasters
						cat("Analysis of randomised branch positions ", s, "\n", sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
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
							}
					}			
				if (branchRandomisation1 == TRUE)
					{
						simRasters = list()
						simRasters = hullRasters
						cat("Analysis of randomised branch positions ", s, "\n", sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
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
							}
					}		
				if (impactOnVelocity == TRUE)
					{
						if (distPermutations == TRUE)
							{
								cat("Analysis on environmental distances permutation ", s, "\n", sep="")
								for (h in 1:length(envVariables))
									{
										for (t in 1:nberOfExtractionFiles)
											{
												distancesSim[[t]][,h] = distances[[t]][sample(dim(distances[[t]])[1]),h]
											}
									}
							}	else		{
								for (t in 1:nberOfExtractionFiles)
									{
										distancesSim[[t]][,1] = distances[[t]][,1]
									}	
								for (h in 2:length(envVariables))
									{
										if (fourCells == TRUE) directions = 4
										if (fourCells == FALSE) directions = 8
										if (resistances[h] == FALSE)
											{
												if (leastCostDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], mean, directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="c", multpl=F, scl=T)
													}
												if (commuteDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], mean, directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
													}
												if (rSPDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], mean, directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
													}
											}	else	{
												if (leastCostDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="c", multpl=F, scl=T)
													}
												if (commuteDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
													}
												if (rSPDistance == TRUE)
													{
														simTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
														simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
													}
											}
										buffer = list()
										# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
										for (t in 1:nberOfExtractionFiles) {
												mat = matrix(nrow=nberOfConnections[t], ncol=1)
												if (straightLineDistance == TRUE)
													{	
														linesList = list()
														for (i in 1:length(fromCoorRand[[t]][,1]))
															{
																points = rbind(fromCoorRand[[t]][i,], toCoorRand[[t]][i,])
																linesList[[i]] = Lines(list(Line(points)),i)
															}
														lines = SpatialLines(linesList)
														extractions = raster::extract(simRasters[[h]], lines)
														for (i in 1:length(fromCoorRand[[t]][,1]))
															{
																mat[i] = sum(extractions[[i]], na.rm=T)
															}	
													}	
												if (leastCostDistance == TRUE)
													{	
														mat[] = diag(costDistance(simTrEnvVariableCorr, fromCoorRand[[t]], toCoorRand[[t]]))
													}
												if (randomWalkDistance == TRUE)
													{
														branchesNotNA = which(!((is.na(raster::extract(simRasters[[h]],fromCoorRand[[t]][])))|
																			   (is.na(raster::extract(simRasters[[h]],toCoorRand[[t]][])))))
														simRasterName = paste("CS_rasters/",names(simRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
														if (juliaCSImplementation == FALSE)
															{
																mat[branchesNotNA,] = circuitScape1(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], fourCells, 
																									fromCoorRand[[t]][branchesNotNA,], toCoorRand[[t]][branchesNotNA,], OS, outputName,
																									t, nberOfCores_CS)
															}	else	{
																mat[branchesNotNA,] = circuitScape2(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], fourCells, 
																									fromCoorRand[[t]][branchesNotNA,], toCoorRand[[t]][branchesNotNA,], OS, outputName,
																									t, nberOfCores_CS)
															}
													}
												if (commuteDistance == TRUE)
													{	
														for (i in 1:length(fromCoorRand[[t]][,1]))
															{
																spatialPoints = SpatialPoints(cbind(c(fromCoorRand[[t]][i,1],toCoorRand[[t]][i,1]),
																									c(fromCoorRand[[t]][i,2],toCoorRand[[t]][i,2])))
																mat[i] = commuteDistance(simTrEnvVariableCorr, spatialPoints)
															}
													}
												if (rSPDistance == TRUE)
													{
														mat[] = diag(rSPDistance(simTrEnvVariableCorr, fromCoorRand[[t]], toCoorRand[[t]],
																				 theta=thetaValue, totalNet="total", method=1))
													}	
												colnames(mat) = names(envVariables[[h]])
												buffer[[t]] = mat
												# mat
											}		
										for (t in 1:length(buffer))
											{
												buffer[[t]][!is.finite(buffer[[t]][])] = NA
												buffer[[t]][buffer[[t]][]==-1] = NA
												distancesSim[[t]][,h] = buffer[[t]][]
											}	
									}		
							}
					}
				if (impactOnDirection == TRUE)
					{
						simMeanExtractions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simRateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simMeanDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						for (h in 2:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										envValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
										for (i in 1:length(ancestralNodes))
											{
												ancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
												envValues = envValues + raster::extract(envVariables[[h]],
															cbind(fromCoorRand[[t]][ancestralBranch,1],fromCoorRand[[t]][ancestralBranch,2]))
											}
										simMeanExtractions[t,h] = envValues + mean(raster::extract(envVariables[[h]], cbind(toCoorRand[[t]][,1],toCoorRand[[t]][,2])), na.rm=T)
										if ((!is.na(meanEnvValues[t,h])) & (!is.na(simMeanExtractions[t,h])))
											{
												if (meanEnvValues[t,h] < simMeanExtractions[t,h]) meanEnvValuesHigher[t,h] = meanEnvValuesHigher[t,h] + 1
												if (meanEnvValues[t,h] > simMeanExtractions[t,h]) meanEnvValuesLower[t,h] = meanEnvValuesLower[t,h] + 1
											}	else	{
												meanEnvValuesHigher[t,h] = NA; meanEnvValuesLower[t,h] = NA
											}
										diffs = raster::extract(envVariables[[h]], fromCoorRand[[t]])-raster::extract(envVariables[[h]], toCoorRand[[t]])
										simRateOfPositiveDifferences[t,h] = sum(diffs[!is.na(diffs)] > 0)/length(diffs[!is.na(diffs)])
										simMeanDifferences[t,h] = mean(diffs, na.rm=T)
										if ((!is.na(rateOfPositiveDifferences[t,h])) & (!is.na(sum(diffs > 0))))
											{
												if (rateOfPositiveDifferences[t,h] < simRateOfPositiveDifferences[t,h])
													{
														rateOfPositiveDifferencesHigher[t,h] = rateOfPositiveDifferencesHigher[t,h] + 1
													}
												if (rateOfPositiveDifferences[t,h] > simRateOfPositiveDifferences[t,h])
													{
														rateOfPositiveDifferencesLower[t,h] = rateOfPositiveDifferencesLower[t,h] + 1
													}
											}	else		{
												rateOfPositiveDifferencesLower[t,h] = NA; rateOfPositiveDifferencesHigher[t,h] = NA
											}
										if ((!is.na(meanDifferences[t,h])) & (!is.na(mean(diffs))))
											{
												if (meanDifferences[t,h] < simMeanDifferences[t,h]) meanDifferencesHigher[t,h] = meanDifferencesHigher[t,h] + 1
												if (meanDifferences[t,h] > simMeanDifferences[t,h]) meanDifferencesLower[t,h] = meanDifferencesLower[t,h] + 1
											}	else		{
												meanDifferencesLower[t,h] = NA; meanDifferencesHigher[t,h] = NA
											}
									}
							}
					}
				
				# 3.2. Linear regressions analyses on randomised datasets:	
				if (impactOnVelocity == TRUE)
					{
						simUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRRsquarePValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						if (GLM == TRUE)
							{
								simMultiGLMcoefficients = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
								simMultiGLMcoefficientPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
								simMultiGLMCAuniqueContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
								simMultiGLMCAcommonContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
							}
						if (distPermutations == TRUE)
							{
								for (t in 1:nberOfExtractionFiles)
									{
										distVariables = paste("dispersalTime[[",t,"]]"," ~ distancesSim[[",t,"]][,1]",sep="")
										LM = lm(as.formula(distVariables))
										simUniLRRsquares1[t,1] = summary(LM)$r.squared
										f = summary(LM)$fstatistic
										if (is.numeric(f))
											{
												p = pf(f[1],f[2],f[3],lower.tail=F)
												attributes(p) = NULL; simUniLRRsquarePValues[t,h] = p
											}	else		{
												simUniLRRsquarePValues[t,h] = NA 
											}
									}
							}
						for (h in 1:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										distVariables = paste("dispersalTime[[",t,"]]"," ~ distancesSim[[",t,"]][,",h,"]",sep="")
										LM = lm(as.formula(distVariables))
										simUniLRRsquares1[t,h] = summary(LM)$r.squared
										f = summary(LM)$fstatistic; p = pf(f[1],f[2],f[3],lower.tail=F)
										attributes(p) = NULL; simUniLRRsquarePValues[t,h] = p
										simUniLRDeltaRsquares1[t,h] = simUniLRRsquares1[t,h]-simUniLRRsquares1[t,1]
										nonNaNdispersalTimes = dispersalTime[[t]][!is.na(distancesSim[[t]][,1])]
										nonNaNdistancesSim1 = distancesSim[[t]][,1][!is.na(distancesSim[[t]][,1])]
										nonNaNdistancesSim2 = distancesSim[[t]][,h][!is.na(distancesSim[[t]][,1])]
										distVariables = paste("nonNaNdistancesSim2 ~ nonNaNdistancesSim1",sep="")
										LM1 = lm(as.formula(distVariables)); residuals_LM1 = LM1$residuals; indices = which(!is.na(residuals_LM1))
										distVariables = paste("nonNaNdispersalTimes[indices] ~ nonNaNdistancesSim2[indices] + residuals_LM1",sep="")
										LM2 = lm(as.formula(distVariables)); simUniLRRsquares2[t,h] = summary(LM2)$r.squared
										simUniLRDeltaRsquares2[t,h] = simUniLRRsquares2[t,h]-simUniLRRsquares1[t,1]
									}
							}
						if (GLM == TRUE)
							{
								if (all == FALSE)
									{	
										for (i in 1:length(envVariables))
											{	
												for (t in 1:nberOfExtractionFiles)
													{
														distVariables = paste("dispersalTime ~ ", names(envVariables[[1]]), sep="")
														matMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
														matMultiGLM[,1] = dispersalTime[[t]]
														matMultiGLM[,2] = distances[[t]][,1]
														matMultiGLMNames = c("dispersalTime")
														matMultiGLMNames = c(matMultiGLMNames, names(envVariables[[1]]))
														for (h in 2:length(envVariables))
															{
																matMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
																distVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
																if (h == i)
																	{
																		matMultiGLM[,h+1] = distancesSim[[t]][,h]
																	}	else	{
																		matMultiGLM[,h+1] = distances[[t]][,h]
																	}
															}
														matMultiGLM = as.data.frame(matMultiGLM)
														names(matMultiGLM) = matMultiGLMNames
														m = min(unlist(lapply(as.matrix(matMultiGLM),min)), na.rm=T)	
														matMultiGLM = lapply(matMultiGLM, preLogTransformation, m)
														matMultiGLM = lapply(matMultiGLM, logTransformation)
														matMultiGLM = lapply(matMultiGLM, zTransformation)
														matMultiGLM = as.data.frame(matMultiGLM)
														names(matMultiGLM) = matMultiGLMNames
														form = as.formula(distVariables)
														multiGLM = stats::glm(form, data=matMultiGLM)
														if (commonalityAnalysis == TRUE)
															{
																CA = calc.yhat(multiGLM, prec=5)
															}
														for (h in 1:length(envVariables))
															{
																names(multiGLM$coefficients)[h] = names(envVariables[[h]])
															}
														simMultiGLMcoefficients[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
														simMultiGLMcoefficientPValues[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
														if (commonalityAnalysis == TRUE)
															{
																simMultiGLMCAuniqueContributions[t,i] = CA$PredictorMetrics[h,"Unique"]
																simMultiGLMCAcommonContributions[t,i] = CA$PredictorMetrics[h,"Common"]
															}
													}
											}
									}	else		{
										for (t in 1:nberOfExtractionFiles)
											{
												distVariables = paste("dispersalTime ~ ", names(envVariables[[1]]), sep="")
												matMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
												matMultiGLM[,1] = dispersalTime[[t]]
												matMultiGLM[,2] = distances[[t]][,1]
												matMultiGLMNames = c("dispersalTime")
												matMultiGLMNames = c(matMultiGLMNames, names(envVariables[[1]]))
												for (h in 2:length(envVariables))
													{
														matMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
														distVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
														# distVariables = paste(distVariables, " + ", "distancesSim[[", t, "]][,", h, "]", sep="")
														matMultiGLM[,h+1] = distancesSim[[t]][,h]	
													}
												matMultiGLM = as.data.frame(matMultiGLM)
												names(matMultiGLM) = matMultiGLMNames
												matMultiGLM = lapply(matMultiGLM, zTransformation)
												matMultiGLM = as.data.frame(matMultiGLM)
												names(matMultiGLM) = matMultiGLMNames
												form = as.formula(distVariables)
												multiGLM = stats::glm(form, data=matMultiGLM)
												if (commonalityAnalysis == TRUE)
													{
														CA = calc.yhat(multiGLM, prec=5)
													}
												for (h in 1:length(envVariables))
													{
														names(multiGLM$coefficients)[h] = names(envVariables[[h]])
													}
												simMultiGLMcoefficients[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
												simMultiGLMcoefficientPValues[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
												if (commonalityAnalysis == TRUE)
													{
														simMultiGLMCAuniqueContributions[t,i] = CA$PredictorMetrics[h,"Unique"]
														simMultiGLMCAcommonContributions[t,i] = CA$PredictorMetrics[h,"Common"]
													}
											}
									}
							}
						if (distPermutations == TRUE)
							{
								H = 1
							}	else		{
								H = 2
							}
					}
				if (impactOnDirection == TRUE) H = 2

				# 3.3.1. Count of lower & higher values compared to real/observed values for this randomisation:			
				if (impactOnVelocity == TRUE)
					{
						for (h in H:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										if (simUniLRRsquares1[t,h] < realUniLRRsquares1[t,h])
											{
												uniLRRsquaresLower1[t,h] = uniLRRsquaresLower1[t,h] + 1;
											}	else	{
												uniLRRsquaresHigher1[t,h] = uniLRRsquaresHigher1[t,h] + 1;
											}
										if (simUniLRRsquarePValues[t,h] < realUniLRRsquarePValues[t,h])
											{
												uniLRRsquarePValuesLower[t,h] = uniLRRsquarePValuesLower[t,h] + 1;
											}	else	{
												uniLRRsquarePValuesHigher[t,h] = uniLRRsquarePValuesHigher[t,h] + 1;
											}
										if (simUniLRDeltaRsquares1[t,h] < realUniDeltaRsquares1[t,h])
											{
												uniLRdeltaRsquaresLower1[t,h] = uniLRdeltaRsquaresLower1[t,h] + 1;
											}	else	{
												uniLRdeltaRsquaresHigher1[t,h] = uniLRdeltaRsquaresHigher1[t,h] + 1;
											}
										if (simUniLRRsquares2[t,h] < realUniLRRsquares2[t,h])
											{
												uniLRRsquaresLower2[t,h] = uniLRRsquaresLower2[t,h] + 1;
											}	else	{
												uniLRRsquaresHigher2[t,h] = uniLRRsquaresHigher2[t,h] + 1;
											}
										if (simUniLRDeltaRsquares2[t,h] < realUniDeltaRsquares2[t,h])
											{
												uniLRdeltaRsquaresLower2[t,h] = uniLRdeltaRsquaresLower2[t,h] + 1;
											}	else	{
												uniLRdeltaRsquaresHigher2[t,h] = uniLRdeltaRsquaresHigher2[t,h] + 1;
											}
										if (GLM == TRUE)
											{
												if (!is.na(simMultiGLMcoefficients[t,h]) & !is.na(realMultiGLMcoefficients[t,h]))
													{	
														if (simMultiGLMcoefficients[t,h] < realMultiGLMcoefficients[t,h])
															{
																multiGLMcoefficientsLower[t,h] = multiGLMcoefficientsLower[t,h] + 1
															}	else	{
																multiGLMcoefficientsHigher[t,h] = multiGLMcoefficientsHigher[t,h] + 1
															}
													}									
												if (!is.na(simMultiGLMcoefficientPValues[t,h]) & !is.na(realMultiGLMcoefficientPValues[t,h]))
													{	
														if (simMultiGLMcoefficientPValues[t,h] < realMultiGLMcoefficientPValues[t,h])
															{
																multiGLMcoefficientPValuesLower[t,h] = multiGLMcoefficientPValuesLower[t,h] + 1;
															}	else	{
																multiGLMcoefficientPValuesHigher[t,h] = multiGLMcoefficientPValuesHigher[t,h] + 1;
															}
													}
												if (commonalityAnalysis == TRUE)
													{	
														if (simMultiGLMCAcommonContributions[t,h] < realMultiGLMCAcommonContributions[t,h])
															{
																multiGLMCAcommonContributionsLower[t,h] = multiGLMCAcommonContributionsLower[t,h] + 1;
															}	else	{
																multiGLMCAcommonContributionsHigher[t,h] = multiGLMCAcommonContributionsHigher[t,h] + 1;
															}
														if (simMultiGLMCAuniqueContributions[t,h] < realMultiGLMCAuniqueContributions[t,h])
															{
																multiGLMCAuniqueContributionsLower[t,h] = multiGLMCAuniqueContributionsLower[t,h] + 1;
															}	else	{
																multiGLMCAuniqueContributionsHigher[t,h] = multiGLMCAuniqueContributionsHigher[t,h] + 1;
															}
													}
											}
									}
							}
					}
					
				# 3.3.2. Computations of the median and BF values for this randomisation:	
				if (impactOnVelocity == TRUE)
					{
						for (h in H:length(envVariables))
							{
								uniLRmedianRsquaresSim[s,h] = median(simUniLRRsquares1[,h])
								uniLRmedianRsquarePValuesSim[s,h] = median(simUniLRRsquarePValues[,h])
								uniLRmedianDeltaRsquaresSim[s,h] = median(simUniLRRsquares1[,h]-simUniLRRsquares1[,1])
								if (GLM == TRUE)
									{
										multiGLMmedianCoefficientsSim[s,h] = median(simMultiGLMcoefficients[,h])
										multiGLMmedianCoefficientPValuesSim[s,h] = median(simMultiGLMcoefficientPValues[,h])
										if (commonalityAnalysis == TRUE)
											{
												multiGLMmedianCAuniqueContributionsSim[s,h] = median(simMultiGLMCAuniqueContributions[,h])
												multiGLMmedianCAcommonContributionsSim[s,h] = median(simMultiGLMCAcommonContributions[,h])
											}
									}
							}
						for (h in H:length(envVariables))
							{
								c = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if (simUniLRDeltaRsquares1[t,h] < realUniDeltaRsquares1[t,h]) c = c+1
									}
								f = c/nberOfExtractionFiles; bf = f/(1-f)
								uniLRdeltaRsquaresRandomisationBFs1[h,s] = round(bf, 4)
								c = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if (simUniLRDeltaRsquares2[t,h] < realUniDeltaRsquares2[t,h]) c = c+1
									}
								f = c/nberOfExtractionFiles; bf = f/(1-f)
								uniLRdeltaRsquaresRandomisationBFs2[h,s] = round(bf, 4)
							}
						if (showingPlots == TRUE) dev.off()
					}
				if (impactOnDirection == TRUE)
					{
						for (h in H:length(envVariables))
							{
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simMeanExtractions[t,h]))&(!is.na(meanEnvValues[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (simMeanExtractions[t,h] > meanEnvValues[t,h]) c = c+1
													}	else		{
														if (simMeanExtractions[t,h] < meanEnvValues[t,h]) c = c+1
													}
											}	else		{
												missingValues = missingValues + 1
											}
									}
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								meanEnvValuesRandomisationBFs[h,s] = round(bf, 4)
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simMeanDifferences[t,h]))&(!is.na(meanDifferences[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (simMeanDifferences[t,h] < meanDifferences[t,h]) c = c+1
													}	else		{
														if (simMeanDifferences[t,h] > meanDifferences[t,h]) c = c+1
													}
											}	else		{
												missingValues = missingValues + 1
											}
									}
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								meanDifferencesRandomisationBFs[h,s] = round(bf, 4)
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simRateOfPositiveDifferences[t,h])) & (!is.na(rateOfPositiveDifferences[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (simRateOfPositiveDifferences[t,h] < rateOfPositiveDifferences[t,h]) c = c+1
													}	else		{
														if (simRateOfPositiveDifferences[t,h] > rateOfPositiveDifferences[t,h]) c = c+1
													}
											}	else		{
												missingValues = missingValues + 1
											}
									}	
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								rateOfPositiveDifferencesRandomisationBFs[h,s] = round(bf, 4)
							}
					}
			} # End of randomisations
				
		# 3.4.1. Computations of p-values for the histogram constructions:	
		if (impactOnVelocity == TRUE)
			{	
				for (h in H:length(envVariables))
					{
						for (t in 1:nberOfExtractionFiles)
							{
								uniLRRsquaresRandomisationPValues1[t,h] = (nberOfRandomisations-uniLRRsquaresLower1[t,h])/nberOfRandomisations
								uniLRRsquarePValuesRandomisationPValues[t,h] = (uniLRRsquarePValuesLower[t,h])/nberOfRandomisations
								uniLRdeltaRsquaresRandomisationPValues1[t,h] = (nberOfRandomisations-uniLRdeltaRsquaresLower1[t,h])/nberOfRandomisations
								uniLRRsquaresRandomisationPValues2[t,h] = (nberOfRandomisations-uniLRRsquaresLower2[t,h])/nberOfRandomisations
								uniLRdeltaRsquaresRandomisationPValues2[t,h] = (nberOfRandomisations-uniLRdeltaRsquaresLower2[t,h])/nberOfRandomisations
								if (GLM == TRUE)
									{
										multiGLMcoefficientsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMcoefficientsLower[t,h])/nberOfRandomisations
										multiGLMcoefficientPValuesRandomisationPValues[t,h] = (multiGLMcoefficientPValuesLower[t,h])/nberOfRandomisations
										if (commonalityAnalysis == TRUE)
											{
									multiGLMCAuniqueContributionsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMCAuniqueContributionsLower[t,h])/nberOfRandomisations
									multiGLMCAcommonContributionsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMCAcommonContributionsLower[t,h])/nberOfRandomisations
											}
									}
							}
					}
			}	
		if (impactOnDirection == TRUE)
			{	
				for (h in H:length(envVariables))
					{
						for (t in 1:nberOfExtractionFiles)
							{
								if (resistances[h] == TRUE)
									{
										meanEnvValuesRandomisationPValues[t,h] = (meanEnvValuesLower[t,h])/nberOfRandomisations
										rateOfPositiveDifferencesRandomisationPValues[t,h] = (rateOfPositiveDifferencesHigher[t,h])/nberOfRandomisations
										meanDifferencesRandomisationPValues[t,h] = (meanDifferencesHigher[t,h])/nberOfRandomisations	
									}	else		{
										meanEnvValuesRandomisationPValues[t,h] = (meanEnvValuesHigher[t,h])/nberOfRandomisations
										rateOfPositiveDifferencesRandomisationPValues[t,h] = (rateOfPositiveDifferencesLower[t,h])/nberOfRandomisations
										meanDifferencesRandomisationPValues[t,h] = (meanDifferencesLower[t,h])/nberOfRandomisations						
									}
							}
					}
			}
		
		# 3.4.2. Computations of the "median" p-values (one-tailed test, except for multiGLM and CA coefficients: two-tailed tests):		
		if (impactOnVelocity == TRUE)
			{
				uniLRmedianRsquaresLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
				uniLRmedianRsquaresPValues = matrix(nrow=1, ncol=length(envVariables))
				uniLRmedianRsquarePValuesLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
				uniLRmedianRsquarePValuesPValues = matrix(nrow=1, ncol=length(envVariables))
				uniLRmedianDeltaRsquaresLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
				uniLRmedianDeltaRsquaresPValues = matrix(nrow=1, ncol=length(envVariables))
				if (GLM == TRUE)
					{
						multiGLMmedianCoefficientsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))	
						multiGLMmedianCoefficientsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
						multiGLMmedianCoefficientsPValues = matrix(nrow=1, ncol=length(envVariables))
						multiGLMmedianCoefficientPValuesLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
						multiGLMmedianCoefficientPValuesPValues = matrix(nrow=1, ncol=length(envVariables))
						multiGLMmedianCAuniqueContributionsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))	
						multiGLMmedianCAuniqueContributionsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
						multiGLMmedianCAuniqueContributionsPValues = matrix(nrow=1, ncol=length(envVariables))
						multiGLMmedianCAcommonContributionsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))	
						multiGLMmedianCAcommonContributionsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
						multiGLMmedianCAcommonContributionsPValues = matrix(nrow=1, ncol=length(envVariables))
					}
				for (h in H:length(envVariables))
					{
						for (s in 1:nberOfRandomisations)
							{
								if (uniLRmedianRsquaresSim[s,h] < uniLRmedianRsquares[h])
									{
										uniLRmedianRsquaresLowerValues[h] = uniLRmedianRsquaresLowerValues[h] + 1
									}
								if (uniLRmedianRsquarePValuesSim[s,h] < uniLRmedianRsquarePValues[h])
									{
										uniLRmedianRsquarePValuesLowerValues[h] = uniLRmedianRsquarePValuesLowerValues[h] + 1
									}
								if (uniLRmedianDeltaRsquaresSim[s,h] < uniLRmedianDeltaRsquares[h])
									{
										uniLRmedianDeltaRsquaresLowerValues[h] = uniLRmedianDeltaRsquaresLowerValues[h] + 1
									}
								if (GLM == TRUE)
									{
										if (!is.na(multiGLMmedianCoefficients[h]))
											{	
												if (multiGLMmedianCoefficientsSim[s,h] < multiGLMmedianCoefficients[h])
													{
														multiGLMmedianCoefficientsLowerValues[h] = multiGLMmedianCoefficientsLowerValues[h] + 1
													}	else	{
														multiGLMmedianCoefficientsHigherValues[h] = multiGLMmedianCoefficientsHigherValues[h] + 1
													}
												if (multiGLMmedianCoefficientPValuesSim[s,h] < multiGLMmedianCoefficientPValues[h])
													{
														multiGLMmedianCoefficientPValuesLowerValues[h] = multiGLMmedianCoefficientPValuesLowerValues[h] + 1
													}
											}
										if (commonalityAnalysis == TRUE)
											{	
												if (multiGLMmedianCAuniqueContributionsSim[s,h] < multiGLMmedianCAuniqueContributions[h])
													{
														multiGLMmedianCAuniqueContributionsLowerValues[h] = multiGLMmedianCAuniqueContributionsLowerValues[h] + 1
													}	else	{
														multiGLMmedianCAuniqueContributionsHigherValues[h] = multiGLMmedianCAuniqueContributionsHigherValues[h] + 1
													}
												if (multiGLMmedianCAcommonContributionsSim[s,h] < multiGLMmedianCAcommonContributions[h])
													{
														multiGLMmedianCAcommonContributionsLowerValues[h] = multiGLMmedianCAcommonContributionsLowerValues[h] + 1
													}	else	{
														multiGLMmedianCAcommonContributionsHigherValues[h] = multiGLMmedianCAcommonContributionsHigherValues[h] + 1
													}
											}
									}				
							}
						uniLRmedianRsquaresPValues[h] = (nberOfRandomisations-uniLRmedianRsquaresLowerValues[h])/nberOfRandomisations
						uniLRmedianRsquarePValuesPValues[h] = (uniLRmedianRsquarePValuesLowerValues[h])/nberOfRandomisations
						uniLRmedianDeltaRsquaresPValues[h] = (nberOfRandomisations-uniLRmedianDeltaRsquaresPValues[h])/nberOfRandomisations
						if (GLM == TRUE)
							{
								multiGLMmedianCoefficientsPValues[h] = (nberOfRandomisations-multiGLMmedianCoefficientsLowerValues[h])/nberOfRandomisations
								multiGLMmedianCoefficientPValuesPValues[h] = (multiGLMmedianCoefficientPValuesLowerValues[h])/nberOfRandomisations
								if (commonalityAnalysis == TRUE)
									{
								multiGLMmedianCAuniqueContributionsPValues[h] = (nberOfRandomisations-multiGLMmedianCAuniqueContributionsLowerValues[h])/nberOfRandomisations
								multiGLMmedianCAcommonContributionsPValues[h] = (nberOfRandomisations-multiGLMmedianCAcommonContributionsLowerValues[h])/nberOfRandomisations
									}
							}
					}
			}
		
		# 3.5. Creation of the pdf file containing histograms:
		if (impactOnVelocity == TRUE)
			{
				fileName = paste(outputName, "_randomisation_results.pdf", sep="")	
				if ((plottingHistograms == TRUE)&(nberOfExtractionFiles > 49)&(nberOfRandomisations > 1))
					{
						if (GLM == TRUE)
							{
								if (commonalityAnalysis == TRUE)
									{
										pdf(fileName, width=(4*length(envVariables)), height=(5*4))
										par(mfrow=c(7,length(envVariables)))
									}	else		{
										pdf(fileName, width=(4*length(envVariables)), height=(3*4))
										par(mfrow=c(5,length(envVariables)))
									}
							}	else	{
								pdf(fileName, width=(4*length(envVariables)), height=(2.25*4))
								par(mfrow=c(3,length(envVariables)))
							}
						breakList = (0:50)/50
						xMin = min(realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,1], na.rm=T)	
						for (h in 2:length(envVariables))
							{
								if (xMin > min(realUniLRRsquares1[,h], na.rm=T)) xMin = min(realUniLRRsquares1[,h], na.rm=T)
								if (xMax < max(realUniLRRsquares1[,h], na.rm=T)) xMax = max(realUniLRRsquares1[,h], na.rm=T)

							}
						for (h in 1:length(envVariables))
							{
								min = min(realUniLRRsquares1[,h]); max = max(realUniLRRsquares1[,h])
								hist(realUniLRRsquares1[,h], freq=T, xlim=c(0,xMax), breaks=seq(0,xMax,by=(xMax-0)/50), xlab="Univariate LR R2's", 
									 main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
							}
						if (distPermutations == FALSE)
							{	
								yMax = 0	
								for (h in 2:length(envVariables))
									{
										a = hist(uniLRRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
										if (yMax < max(a$counts)) yMax = max(a$counts)
									}
								plot.new()
							}	else	{
								yMax = 0	
								for (h in 1:length(envVariables))
									{
										a = hist(uniLRRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
										if (yMax < max(a$counts)) yMax = max(a$counts)
									}
								hist(uniLRRsquaresRandomisationPValues1[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
									 xlab="Univariate LR R2 p-values (randomisations)", main="geographical distance", cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)	
							}		
						for (h in 2:length(envVariables))
							{
								hist(uniLRRsquaresRandomisationPValues1[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
									 xlab="Univariate LR R2 p-values (randomisations)", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
							}
						yMax = 0	
						for (h in 2:length(envVariables))
							{
								a = hist(uniLRdeltaRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
								if (yMax < max(a$counts)) yMax = max(a$counts)
							}
						plot.new()
						for (h in 2:length(envVariables))
							{
								hist(uniLRdeltaRsquaresRandomisationPValues1[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
									 xlab="Univariate LR delta R2 (Q) p-values (randomisations)", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
							}
						if (GLM == TRUE)
							{
								xMin = min(realMultiGLMcoefficients[,1], na.rm=T); xMax = max(realMultiGLMcoefficients[,1], na.rm=T)
								for (h in 2:length(envVariables))
									{
										if (xMin > min(realMultiGLMcoefficients[,h], na.rm=T))
											{
												xMin = min(realMultiGLMcoefficients[,h], na.rm=T)
											}
										if (xMax < max(realMultiGLMcoefficients[,h], na.rm=T))
											{
												xMax = max(realMultiGLMcoefficients[,h], na.rm=T)
											}
									}
								for (h in 1:length(envVariables))
									{
										if (is.na(mean(realMultiGLMcoefficients[,h])))
											{	
												plot.new()
											}	else	{
												hist(realMultiGLMcoefficients[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
													 xlab="Multivariate GLM coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}	
									}
								if (distPermutations == FALSE)
									{	
										yMax = 0	
										for (h in 2:length(envVariables))
											{
												a = hist(multiGLMcoefficientsRandomisationPValues[,h], breaks=breakList, plot=F)
												if (yMax < max(a$counts))
													{
														yMax = max(a$counts)
													} 
											}
										plot.new()
									}	else	{
										yMax = 0	
										for (h in 1:length(envVariables))
											{
												a = hist(multiGLMcoefficientsRandomisationPValues[,h], breaks=breakList, plot=F)
												if (yMax < max(a$counts))
													{
														yMax = max(a$counts)
													} 
											}
										hist(multiGLMcoefficientsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
											 xlab="Multivariate GLM coefficient p-values (randomisations)", main="geographical distance",
											 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
									}		
								for (h in 2:length(envVariables))
									{
										hist(multiGLMcoefficientsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
											 xlab="Multivariate GLM coefficient p-values (randomisations)", main=names(envVariables[[h]]),
											 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
									}
								if (commonalityAnalysis == TRUE)
									{	
										if (distPermutations == FALSE)
											{	
												yMax = 0	
												for (h in 2:length(envVariables))
													{
														a = hist(multiGLMCAuniqueContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
														if (yMax < max(a$counts))
															{
																yMax = max(a$counts)
															} 
													}
												plot.new()
											}	else	{
												yMax = 0	
												for (h in 1:length(envVariables))
													{
														a = hist(multiGLMCAuniqueContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
														if (yMax < max(a$counts)) yMax = max(a$counts)
													}
												hist(multiGLMCAuniqueContributionsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
													 xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main="geographical distance", 
													 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}		
										for (h in 2:length(envVariables))
											{
												hist(multiGLMCAuniqueContributionsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
													 xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main=names(envVariables[[h]]), 
													 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}	
										if (distPermutations == FALSE)
											{	
												yMax = 0	
												for (h in 2:length(envVariables))
													{
														a = hist(multiGLMCAcommonContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
														if (yMax < max(a$counts)) yMax = max(a$counts)
													}
												plot.new()
											}	else	{
												yMax = 0	
												for (h in 1:length(envVariables))
													{
														a = hist(multiGLMCAcommonContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
														if (yMax < max(a$counts)) yMax = max(a$counts)
													}
												hist(multiGLMCAcommonContributionsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
													 xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main="geographical distance", 
													 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}		
										for (h in 2:length(envVariables))
											{
												hist(multiGLMCAcommonContributionsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
													 xlab="Multivariate GLM CA common contributions p-values (randomisations)", main=names(envVariables[[h]]), 
													 cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
											}
									}
							}
						dev.off()
					}
			}				
			
		# 3.6. Creation of the text files containing randomisation results:
		if ((nberOfRandomisations > 1) & (impactOnVelocity == TRUE))
			{
				fileName = paste(outputName, "_randomisation_results.txt", sep="")
				if (GLM == TRUE)
					{
						mat = cbind(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
									uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]],
									multiGLMcoefficientsRandomisationPValues[,2:dim(multiGLMcoefficientsRandomisationPValues)[2]])
						if (nberOfExtractionFiles == 1)
							{
								mat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
										  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]], 
										  multiGLMcoefficientsRandomisationPValues[,2:dim(multiGLMcoefficientsRandomisationPValues)[2]]))
							}
						if (commonalityAnalysis == TRUE)
							{
								mat = cbind(mat, multiGLMCAuniqueContributionsRandomisationPValues[,2:dim(multiGLMCAuniqueContributionsRandomisationPValues)[2]],
												 multiGLMCAcommonContributionsRandomisationPValues[,2:dim(multiGLMCAcommonContributionsRandomisationPValues)[2]])
							}
					}	else	{
						mat = cbind(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
									uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]])
						if (alternativeQstat == TRUE)
							{
								mat = cbind(mat, uniLRRsquaresRandomisationPValues2[,2:dim(uniLRRsquaresRandomisationPValues2)[2]], 
											uniLRdeltaRsquaresRandomisationPValues2[,2:dim(uniLRdeltaRsquaresRandomisationPValues2)[2]])
							}
						if (nberOfExtractionFiles == 1)
							{
								mat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
										  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]]))
								if (alternativeQstat == TRUE)
									{
										mat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
										  		  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]],
										  		  uniLRRsquaresRandomisationPValues2[,2:dim(uniLRRsquaresRandomisationPValues2)[2]], 
												  uniLRdeltaRsquaresRandomisationPValues2[,2:dim(uniLRdeltaRsquaresRandomisationPValues2)[2]]))
									}
							}
					}
				uniLRRsquarePValuesNames1 = c(); uniLRdeltaRsquarePValuesNames1 = c(); uniLRRsquarePValuePValuesNames = c()
				uniLRRsquarePValuesNames2 = c(); uniLRdeltaRsquarePValuesNames2 = c()
				if (GLM == TRUE)
					{
						multiGLMcoefficientPValuesNames = c(); multiGLMRsquarePValuePValuesNames = c()
						multiGLMCAuniqueContributionsPValuesNames = c(); multiGLMCAcommonContributionsPValuesNames = c()
					}
				for (h in 2:length(envVariables))
					{
						if (alternativeQstat == TRUE)
							{
								uniLRRsquarePValuesNames1 = cbind(uniLRRsquarePValuesNames1, paste("Uni_LR1_R2_p-values_",names(envVariables[[h]]),sep=""))
								uniLRdeltaRsquarePValuesNames1 = cbind(uniLRdeltaRsquarePValuesNames1, paste("Uni_LR1_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
								uniLRRsquarePValuesNames2 = cbind(uniLRRsquarePValuesNames2, paste("Uni_LR2_R2_p-values_",names(envVariables[[h]]),sep=""))
								uniLRdeltaRsquarePValuesNames2 = cbind(uniLRdeltaRsquarePValuesNames2, paste("Uni_LR2_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
							}	else		{
								uniLRRsquarePValuesNames1 = cbind(uniLRRsquarePValuesNames1, paste("Uni_LR_R2_p-values_",names(envVariables[[h]]),sep=""))
								uniLRdeltaRsquarePValuesNames1 = cbind(uniLRdeltaRsquarePValuesNames1, paste("Uni_LR_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
							}
						if (GLM == TRUE)
							{
								multiGLMcoefficientPValuesNames = cbind(multiGLMcoefficientPValuesNames,
																		paste("Multivariate_GLM_coefficient_p-values_", names(envVariables[[h]]), sep=""))
								if (commonalityAnalysis == TRUE)
									{
										multiGLMCAuniqueContributionsPValuesNames = cbind(multiGLMCAuniqueContributionsPValuesNames, 
										paste("Multivariate_GLM_CA_unique_contributions_", names(envVariables[[h]]), sep=""))
										multiGLMCAcommonContributionsPValuesNames = cbind(multiGLMCAcommonContributionsPValuesNames,
										paste("Multivariate_GLM_CA_common_contributions_", names(envVariables[[h]]), sep=""))
									}
							}
					}
				if (GLM == TRUE)
					{
						names = cbind(uniLRRsquarePValuesNames, uniLRdeltaRsquarePValuesNames, uniLRRsquarePValuePValuesNames, 
									  multiGLMcoefficientPValuesNames, multiGLMRsquarePValuePValuesNames)
						if (commonalityAnalysis == TRUE)
							{
								names = cbind(names, multiGLMCAuniqueContributionsPValuesNames, multiGLMCAcommonContributionsPValuesNames)
							}
					}	else	{
						names = cbind(uniLRRsquarePValuesNames1, uniLRdeltaRsquarePValuesNames1)
						if (alternativeQstat == TRUE) names = cbind(names, uniLRRsquarePValuesNames2, uniLRdeltaRsquarePValuesNames2)
					}
				colnames(mat) = names; write.table(mat, file=fileName, row.names=F, quote=F, sep="\t")
			}
		if ((nberOfRandomisations > 1) & (impactOnDirection == TRUE))
			{
				fileName = paste(outputName, "_randomisation_results.txt", sep="")
				mat = cbind(meanEnvValuesRandomisationPValues[,2:dim(meanEnvValuesRandomisationPValues)[2]], 
							rateOfPositiveDifferencesRandomisationPValues[,2:dim(rateOfPositiveDifferencesRandomisationPValues)[2]])
				if (nberOfExtractionFiles == 1)
					{
						mat = t(c(meanEnvValuesRandomisationPValues[,2:dim(meanEnvValuesRandomisationPValues)[2]], 
								  rateOfPositiveDifferencesRandomisationPValues[,2:dim(rateOfPositiveDifferencesRandomisationPValues)[2]]))
					}
				meanEnvValuesPValuesNames = c(); rateOfPositiveDifferencesPValuesNames = c()
				for (h in 2:length(envVariables))
					{
						meanEnvValuesPValuesNames = cbind(meanEnvValuesPValuesNames, paste("Direction_E_p-values_", names(envVariables[[h]]), sep=""))
						rateOfPositiveDifferencesPValuesNames = cbind(rateOfPositiveDifferencesPValuesNames, paste("Direction_R_p-values_", names(envVariables[[h]]), sep=""))
					}
				names = cbind(meanEnvValuesPValuesNames, rateOfPositiveDifferencesPValuesNames)
				colnames(mat) = names; write.table(mat, file=fileName, row.names=F, quote=F, sep="\t")
			}
		if (nberOfExtractionFiles > 49)
			{
				envVariablesNames = c()
				for (h in 1:length(envVariables))
					{
						envVariablesNames = c(envVariablesNames, names(envVariables[[h]]))
					}
				colNames1 = c(); colNames2 = c()
				for (s in 1:nberOfRandomisations)
					{
						if (alternativeQstat == TRUE)
							{
								colNames1 = c(colNames1, paste("BFs_deltaR2_LR1_randomisation_",s,sep=""))
								colNames2 = c(colNames2, paste("BFs_deltaR2_LR2_randomisation_",s,sep=""))
							}	else	{
								colNames1 = c(colNames1, paste("BFs_deltaR2_randomisation_",s,sep=""))
							}
					}
				if (impactOnVelocity == TRUE)
					{
						row.names(uniLRdeltaRsquaresRandomisationBFs1) = envVariablesNames
						colnames(uniLRdeltaRsquaresRandomisationBFs1) = colNames1
						row.names(uniLRdeltaRsquaresRandomisationBFs2) = envVariablesNames
						colnames(uniLRdeltaRsquaresRandomisationBFs2) = colNames2
						fileName = paste(outputName, "_randomisation_Bayes_factors.txt",sep="")
						if (alternativeQstat == TRUE)
							{
								buffer = cbind(uniLRdeltaRsquaresRandomisationBFs1,uniLRdeltaRsquaresRandomisationBFs2)
							}	else		{
								buffer = uniLRdeltaRsquaresRandomisationBFs1
							}
						write.table(buffer, fileName, quote=F,sep="\t")
					}
				if (impactOnDirection == TRUE)
					{
						if (pathModel == -1)
							{
								row.names(rateOfPositiveDifferencesRandomisationBFs) = envVariablesNames
								colnames(rateOfPositiveDifferencesRandomisationBFs) = "BF"
								fileName = paste(outputName,"_positiveDiffs_Bayes_factors.txt",sep="")
								write.table(rateOfPositiveDifferencesRandomisationBFs, fileName, quote=F, sep="\t")
								row.names(meanDifferencesRandomisationBFs) = envVariablesNames
								colnames(meanDifferencesRandomisationBFs) = "BF"
								fileName = paste(outputName,"_meanDifferences_Bayes_factors.txt",sep="")
								write.table(meanDifferencesRandomisationBFs, fileName, quote=F, sep="\t")
							}
						if (pathModel == 0)
							{
								row.names(meanEnvValuesRandomisationBFs) = envVariablesNames
								colnames(meanEnvValuesRandomisationBFs) = "BF"
								fileName = paste(outputName,"_direction_E_Bayes_factors.txt",sep="")
								write.table(meanEnvValuesRandomisationBFs, fileName, quote=F, sep="\t")
								row.names(rateOfPositiveDifferencesRandomisationBFs) = envVariablesNames
								colnames(rateOfPositiveDifferencesRandomisationBFs) = "BF"
								fileName = paste(outputName,"_direction_R_Bayes_factors.txt",sep="")
								write.table(rateOfPositiveDifferencesRandomisationBFs, fileName, quote=F, sep="\t")
							}
					}
			}
	}
}
