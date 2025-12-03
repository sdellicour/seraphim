spreadFactors = function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), pathModel=2, resistances=list(), avgResistances=list(), 
						 fourCells=FALSE, nberOfRandomisations=0, randomProcedure=3, outputName="", showingPlots=FALSE, nberOfCores=1, 
						 OS="Unix", juliaCSImplementation=FALSE, simulations=FALSE, randomisations=FALSE, minimumConvexHull=TRUE) {

# nberOfCores = 1; OS = "Unix"; juliaCSImplementation = FALSE; simulations = FALSE; randomisations = FALSE; minimumConvexHull = TRUE
impactOnDiffusionVelocity = FALSE; impactOnDispersalLocation = FALSE; impactOnDispersalPath = FALSE
reportingBothQstats = TRUE; savingRandomisationFiles = FALSE
if (pathModel > 0)
	{
		impactOnDiffusionVelocity = TRUE; impactOnDispersalLocation = FALSE; impactOnDispersalPath = FALSE
	}
if (pathModel == 0)
	{
		impactOnDiffusionVelocity = FALSE; impactOnDispersalLocation = TRUE; impactOnDispersalPath = FALSE
	}
if (pathModel < 0)
	{
		impactOnDiffusionVelocity = FALSE; impactOnDispersalLocation = FALSE; impactOnDispersalPath = TRUE
	}
registerDoMC(cores=nberOfCores); nberOfCores_CS = 1
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
nullRaster = envVariables[[1]]; nullRaster[!is.na(nullRaster)] = 1; names(nullRaster) = "null_raster"
newEnvVariables = list(nullRaster); newResistances = c(TRUE); newAvgResistances = c(TRUE)
for (h in 1:length(envVariables))
	{
		newEnvVariables[[h+1]] = envVariables[[h]]
		newEnvVariables[[h+1]][newEnvVariables[[h+1]]<0] = NA
		if(resistances[[h]] == TRUE)
			{
				fric = "_R"
			}	else	{
				fric = "_C"
			}
		names(newEnvVariables[[h+1]]) = paste(names(newEnvVariables[[h+1]]),fric,sep="")
		if (length(resistances) > 0)
			{
				newResistances = c(newResistances, resistances[[h]])
				newAvgResistances = c(newAvgResistances, avgResistances[[h]])
			}
	}
envVariables = newEnvVariables; resistances = newResistances; avgResistances = newAvgResistances
straightLineDistance = FALSE; leastCostDistance = FALSE; circuitscapeDistance = FALSE
externalRandomisations = FALSE; externalSimulations = FALSE; branchRandomisation3 = FALSE
branchRandomisation2 = FALSE; branchRandomisation1 = FALSE
torusRandomisations = FALSE; rastersSimulations = FALSE
if (abs(pathModel) == 1) straightLineDistance = TRUE
if (abs(pathModel) == 2) leastCostDistance = TRUE
if (abs(pathModel) == 3) circuitscapeDistance = TRUE
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
nberOfConnections = rep(NA, nberOfExtractionFiles); totalNberOfConnections = 0
node1 = list(); node2 = list(); fromCoor = list(); toCoor = list()
dispersalTime = list(); treeIDs = list(); datas = list()
if (impactOnDiffusionVelocity == TRUE)
	{
		distances = list()
	}
if (impactOnDispersalLocation == TRUE)
	{
		meanEnvValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		meanEnvDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		rateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		# rateOfPositiveDifferences corresponds to the previous R statistic presented in 2019
	}
for (t in 1:nberOfExtractionFiles)
	{
		if (nchar(localTreesDirectory) == 0)
			{
				fileName = paste(extractionFileName,"_",t,".csv",sep="")
			}	else	{
				fileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
			}	
		data = read.csv(fileName, head=T); datas[[t]] = data; nberOfConnections[t] = dim(data)[1]
		node1[[t]] = matrix(nrow=dim(data)[1], ncol=1); node1[[t]][] = data[,"node1"]
		node2[[t]] = matrix(nrow=dim(data)[1], ncol=1); node2[[t]][] = data[,"node2"]
		fromCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); fromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
		toCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); toCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
		dispersalTime[[t]] = matrix(nrow=dim(data)[1], ncol=1); colnames(dispersalTime[[t]]) = "dispersalTime"
		dispersalTime[[t]][] = data[,"endYear"]-data[,"startYear"]
		totalNberOfConnections = totalNberOfConnections + dim(data)[1]
		if (impactOnDiffusionVelocity == TRUE)
			{
				distances[[t]] = matrix(nrow=dim(data)[1], ncol=length(envVariables))
			}
		if (("treeID"%in%colnames(data)) == TRUE)
			{
				treeIDs[[t]] = data[1,"treeID"]
			}	else	{
				treeIDs[[t]] = "noTreeID"
			}	
	}
hullRasters = list(); hullRasters[1:length(envVariables)] = envVariables[1:length(envVariables)]
if (minimumConvexHull == TRUE)
	{
		points = matrix(nrow=(totalNberOfConnections*2), ncol=2); a = 0
		for (t in 1:nberOfExtractionFiles)
			{
				if (t > 1)
					{
						a = a + nberOfConnections[t-1]
					}
				for (i in 1:nberOfConnections[t])
					{
						index = (a*2) + ((i-1)*2) + 1
						points[index,1] = fromCoor[[t]][i,1]; points[index,2] = fromCoor[[t]][i,2]
						points[(index+1),1] = toCoor[[t]][i,1]; points[(index+1),2] = toCoor[[t]][i,2]								
					}
			}
		points = points[points[,1]>extent(hullRasters[[1]])@xmin,]; points = points[points[,1]<extent(hullRasters[[1]])@xmax,]
		points = points[points[,2]>extent(hullRasters[[1]])@ymin,]; points = points[points[,2]<extent(hullRasters[[1]])@ymax,]
		hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,]); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pointsRaster = rasterize(points, crop(hullRasters[[1]], sps, snap="out")); pointsRaster[!is.na(pointsRaster[])] = 0
	}
for (h in 1:length(envVariables))
	{
		if (minimumConvexHull == TRUE)
			{
				hullRasters[[h]] = crop(hullRasters[[h]], sps, snap="out")
				bufferRaster = hullRasters[[h]]
				hullRasters[[h]] = raster::mask(hullRasters[[h]], sps, snap="out")
				hullRasters[[h]][!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
			}
		names(hullRasters[[h]]) = gsub(".asc","",names(envVariables[[h]]))
		names(hullRasters[[h]]) = gsub(".tif","",names(envVariables[[h]]))
		names(hullRasters[[h]]) = gsub(".gri","",names(envVariables[[h]]))
	}
if (circuitscapeDistance == TRUE)
	{
		extensions = rep("", length(envVariables))
		if ("CS_rasters"%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), "CS_rasters"))
		for (h in 1:length(envVariables))
			{
				extensions[h] = ".asc"
				if (round(res(hullRasters[[h]])[1],10) != round(res(hullRasters[[h]])[2],10)) extensions[h] = ".tif"
				fileName = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
				writeRaster(hullRasters[[h]], fileName, overwrite=T)
			}
	}

# 1. Computation of environmental distances or extraction of environmental values

for (h in 1:length(envVariables))
	{
		if (impactOnDiffusionVelocity == TRUE)
			{
				if (straightLineDistance == TRUE) cat("Computing environmental distances (straight-line path model) for ",names(envVariables[[h]])[1],"\n",sep="")
				if (leastCostDistance == TRUE) cat("Computing environmental distances (least-cost path model) for ",names(envVariables[[h]])[1],"\n",sep="")
				if (circuitscapeDistance == TRUE) cat("Computing environmental distances (Circuitscape path model) for ",names(envVariables[[h]])[1],"\n",sep="")	
				if (fourCells == TRUE) directions = 4
				if (fourCells == FALSE) directions = 8
				if (leastCostDistance == TRUE)
					{
						if (resistances[h] == FALSE)
							{
								trEnvVariable = transition(hullRasters[[h]], mean, directions)
								trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
							}
						if (resistances[h] == TRUE)
							{
								trEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
								trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
							}
					}
				buffer = list()
				buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
				# for (t in 1:nberOfExtractionFiles) {
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
						if (circuitscapeDistance == TRUE)
							{
								envVariableName = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
								branchesNotNA = which(!((is.na(raster::extract(hullRasters[[h]], fromCoor[[t]][])))
													  |(is.na(raster::extract(hullRasters[[h]], toCoor[[t]][])))))
								fromCoor_temp = fromCoor[[t]][branchesNotNA,]; toCoor_temp = toCoor[[t]][branchesNotNA,]
								if (juliaCSImplementation == FALSE)
									{
										mat[branchesNotNA,1] = circuitScape1(hullRasters[[h]], envVariableName, resistances[[h]], avgResistances[[h]],
																			 fourCells, fromCoor_temp, toCoor_temp, OS, outputName, t, nberOfCores_CS)
									}	else	{
										mat[branchesNotNA,1] = circuitScape2(hullRasters[[h]], envVariableName, resistances[[h]], avgResistances[[h]],
																			 fourCells, fromCoor_temp, toCoor_temp, OS, outputName, t, nberOfCores_CS)
									}
							}			
						colnames(mat) = names(envVariables[[h]])
						# buffer[[t]] = mat
						mat
					}		
				for (t in 1:length(buffer))
					{
						buffer[[t]][!is.finite(buffer[[t]][])] = NA # least-cost case
						buffer[[t]][buffer[[t]][]==-1] = NA # CircuitScape case
						buffer[[t]][buffer[[t]][]==-777] = NA # NA value in Circuitscape
						distances[[t]][,h] = buffer[[t]][]
					}
			}
		if (impactOnDispersalLocation == TRUE)
			{
				for (t in 1:nberOfExtractionFiles)
					{
						envValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
						for (i in 1:length(ancestralNodes))
							{
								ancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
								envValues = envValues + raster::extract(envVariables[[h]], cbind(fromCoor[[t]][ancestralBranch,1],fromCoor[[t]][ancestralBranch,2]))
							}
						meanEnvValues[t,h] = mean(c(envValues,raster::extract(envVariables[[h]], toCoor[[t]])), na.rm=T)
						diffs1 = raster::extract(envVariables[[h]], fromCoor[[t]])-raster::extract(envVariables[[h]], toCoor[[t]])
						rateOfPositiveDifferences[t,h] = sum(diffs1[!is.na(diffs1)] > 0)/length(diffs1[!is.na(diffs1)])
						diffs2 = raster::extract(envVariables[[h]], toCoor[[t]])-raster::extract(envVariables[[h]], fromCoor[[t]])
						meanEnvDifferences[t,h] = mean(diffs2, na.rm=T)
					}
				colNames = c()
				for (h in 1:length(envVariables)) colNames = cbind(colNames, names(envVariables[[h]]))
				buffer = meanEnvValues; colnames(buffer) = colNames
				write.csv(buffer, paste(outputName,"_mean_evironmental_values_at_node_locations.csv",sep=""), row.names=F, quote=F)
			}
	}
envVariableNames = names(envVariables[[1]])
for (h in 2:length(envVariables))
	{
		envVariableNames = cbind(envVariableNames, names(envVariables[[h]]))
	}
if (impactOnDiffusionVelocity == TRUE)
	{
		for (t in 1:nberOfExtractionFiles)
			{
				colnames(distances[[t]]) = envVariableNames
			}
	}
if ((nberOfExtractionFiles == 1)&(impactOnDiffusionVelocity == TRUE))
	{
		mat = cbind(dispersalTime[[1]][], distances[[1]][,1:length(envVariables)])
		columnNames = "dispersal_times"
		for (h in 1:length(envVariables))
			{
				columnNames = cbind(columnNames, names(envVariables[[h]]))
			}
		colnames(mat) = columnNames
		write.table(mat, file=paste(outputName,"_env_distances.txt",sep=""), row.names=F, quote=F, sep="\t")
	}
if ((file.exists(outputName))&(impactOnDiffusionVelocity == TRUE))
	{	
		for (t in 1:nberOfExtractionFiles)
			{
				mat = cbind(dispersalTime[[t]][], distances[[t]][,1:length(envVariables)])
				columnNames = cbind("dispersal_times")
				for (h in 1:length(envVariables))
					{
						columnNames = cbind(columnNames, names(envVariables[[h]]))
					}
				colnames(mat) = columnNames
				write.table(mat, file=paste(outputName,"/",outputName, "_tree",t,"_env_distances.txt",sep=""), row.names=F, quote=F, sep="\t")
			}
	}

# 2. Linear regressions of dispersal durations vs environmental distances

if (impactOnDiffusionVelocity == TRUE)
	{
		realUniLRcoefficients1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRcoefficients2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		colNames = c(); for (t in 1:nberOfExtractionFiles) colNames = c(colNames, paste("tree_",treeIDs[[t]],sep=""))
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
				for (t in 1:nberOfExtractionFiles)
					{
						durationTimes = dispersalTime[[t]]; distances1 = distances[[t]][,h]
						LM1 = lm(as.formula("durationTimes ~ distances1"))
						realUniLRcoefficients1[t,h] = summary(LM1)$coefficients[2,"Estimate"]
						realUniLRRsquares1[t,h] = summary(LM1)$r.squared
						durationTimes4 = 4*dispersalTime[[t]]; squareDistances = distances[[t]][,h]^2
						LM2 = lm(as.formula("durationTimes4 ~ squareDistances"))
						realUniLRcoefficients2[t,h] = summary(LM2)$coefficients[2,"Estimate"]
						realUniLRRsquares2[t,h] = summary(LM2)$r.squared
						# plot(distances1, durationTimes); plot(squareDistances, durationTimes4)
					}
			}
		for (h in 2:length(envVariables))
			{
				realUniDeltaRsquares1[,h] = (realUniLRRsquares1[,h]-realUniLRRsquares1[,1])
				realUniDeltaRsquares2[,h] = (realUniLRRsquares2[,h]-realUniLRRsquares2[,1])
			}
		if (nberOfExtractionFiles > 1)
			{
				uniLRcoefficientsNames1 = c(); uniLRRsquaresNames1 = c(); uniLRdeltaRsquaresNames1 = c(); uniLRcoefficientsNames2 = c(); uniLRRsquaresNames2 = c(); uniLRdeltaRsquaresNames2 = c()
				mat = cbind(realUniLRcoefficients1[,1:dim(realUniLRcoefficients1)[2]], realUniLRRsquares1[,1:dim(realUniLRRsquares1)[2]], realUniDeltaRsquares1[,2:dim(realUniDeltaRsquares1)[2]])
				if (reportingBothQstats == TRUE)
					{
						mat = cbind(mat, realUniLRcoefficients2[,1:dim(realUniLRcoefficients2)[2]], realUniLRRsquares2[,1:dim(realUniLRRsquares2)[2]], realUniDeltaRsquares2[,2:dim(realUniDeltaRsquares2)[2]])
					}
				for (h in 1:length(envVariables))
					{
						if (reportingBothQstats == TRUE)
							{
								uniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("LR1_coefficients_",names(envVariables[[h]]),sep=""))	
								uniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("LR1_R2_",names(envVariables[[h]]),sep=""))				
								uniLRcoefficientsNames2 = cbind(uniLRcoefficientsNames2, paste("LR2_coefficients_",names(envVariables[[h]]),sep=""))	
								uniLRRsquaresNames2 = cbind(uniLRRsquaresNames2, paste("LR2_R2_",names(envVariables[[h]]),sep=""))	
								if (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1, paste("LR1_Q_",names(envVariables[[h]]),sep=""))
								if (h > 1) uniLRdeltaRsquaresNames2 = cbind(uniLRdeltaRsquaresNames2, paste("LR2_Q_",names(envVariables[[h]]),sep=""))
							}	else	{
								uniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("LR_coefficients_",names(envVariables[[h]]),sep=""))	
								uniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("LR_R2_",names(envVariables[[h]]), sep=""))				
								if (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1,paste("LR_Q_",names(envVariables[[h]]),sep=""))
							}
					}
				if (reportingBothQstats == TRUE)
					{
						names = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1, uniLRcoefficientsNames2, uniLRRsquaresNames2, uniLRdeltaRsquaresNames2)
					}	else	{
						names = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1)
					}
				colnames(mat) = names; write.table(mat, file=paste(outputName,"_linear_regression_results.txt",sep=""), row.names=F, quote=F, sep="\t")
			}
	}

# 3. Randomisation step (randomisation of phylogenetic branch positions)

if (nberOfRandomisations > 0)
	{	
		if (impactOnDiffusionVelocity == TRUE)
			{
				uniLRdeltaRsquaresRandomisationBFs1 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)	
				uniLRdeltaRsquaresRandomisationBFs2 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)	
			}
		if (impactOnDispersalLocation == TRUE)
			{
				meanEnvValuesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
				meanEnvDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
				rateOfPositiveDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
			}
		for (s in 1:nberOfRandomisations)
			{
				if (impactOnDiffusionVelocity == TRUE)
					{
						uniLRdeltaRsquaresRand = matrix(0, nrow=nberOfExtractionFiles, ncol=nberOfRandomisations)
						distancesSim = list()
						for (t in 1:nberOfExtractionFiles)
							{
								distancesSim[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
							}
					}
				
				# 3.1. Randomisation step
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
						nberOfConnections = rep(NA, nberOfExtractionFiles)
						node1 = list(); node2 = list(); fromCoor = list(); toCoor = list()
						dispersalTime = list(); treeIDs = list()
						for (t in 1:nberOfExtractionFiles)
							{
								if (nchar(localTreesDirectory) == 0)
									{
										fileName = paste(extractionFileName,"_",t,".csv",sep="")
									}	else	{
										fileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv",sep="")
									}	
								data = read.csv(fileName, head=T); nberOfConnections[t] = dim(data)[1]		
								node1[[t]] = matrix(nrow=dim(data)[1], ncol=1); node1[[t]][] = data[,"node1"]
								node2[[t]] = matrix(nrow=dim(data)[1], ncol=1); node2[[t]][] = data[,"node2"]
								fromCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); fromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
								toCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); toCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
								dispersalTime[[t]] = matrix(nrow=dim(data)[1], ncol=1); colnames(dispersalTime[[t]]) = "dispersalTime"
								dispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
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
						cat("Analysis on raster simulation ",s,"\n",sep="")
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
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
							{
								counter1 = 0; twoPointsOnTheGrid = FALSE
								while (twoPointsOnTheGrid == FALSE)
									{
										twoPointsOnTheGrid = TRUE; counter1 = counter1+1
										fromCoorRand[[t]][,] = NA; toCoorRand[[t]][,] = NA
										if (showingPlots == TRUE)
											{
												if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
												if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
												lines(points[hull,], lwd=0.5, col="black")
												text1 = paste("randomisation of branch positions, sampled tree ",t,sep="")
												mtext(text1, col="black", cex=0.7, line=0)
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
																										pt2_rotated = pt1
																										pt2NAarea = FALSE
																										twoPointsOnTheGrid = FALSE
																									}	else	{
																										pt2_rotated = pt1
																										pt2NAarea = FALSE
																										twoPointsOnTheGrid = TRUE
																									}
																							}
																						if (ancestralNodeID == TRUE)
																							{
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
								if (savingRandomisationFiles)
									{
										temp = datas[[t]]; temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
										write.csv(temp, paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)						
									}
							}
					}
				if (branchRandomisation2 == TRUE)
					{
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
							{
								if (showingPlots == TRUE)
									{
										if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
										if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
										lines(points[hull,], lwd=0.5, col="black")
										text1 = paste("randomisation of branch positions, sampled tree ",t,sep="")
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
								if (savingRandomisationFiles)
									{
										temp = datas[[t]]; temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
										write.csv(temp, paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)						
									}
							}
					}			
				if (branchRandomisation1 == TRUE)
					{
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor
						for (t in 1:nberOfExtractionFiles)
							{
								if (showingPlots == TRUE)
									{
										if (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
										if (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=T)
										plot(sps, lwd=0.5, border="black", add=T)
										text1 = paste("randomisation of branch positions, sampled tree ",t,sep="")
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
								if (savingRandomisationFiles)
									{
										temp = datas[[t]]; temp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
										write.csv(temp, paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)						
									}
							}
					}		
				if (impactOnDiffusionVelocity == TRUE)
					{
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
									}	else	{
										if (leastCostDistance == TRUE)
											{
												simTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
												simTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="c", multpl=F, scl=T)
											}
									}
								buffer = list()
								buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
								# for (t in 1:nberOfExtractionFiles) {
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
										if (circuitscapeDistance == TRUE)
											{
												simRasterName = paste("CS_rasters/",names(simRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
												branchesNotNA = which(!((is.na(raster::extract(simRasters[[h]],fromCoorRand[[t]][])))|
																	   (is.na(raster::extract(simRasters[[h]],toCoorRand[[t]][])))))
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
										colnames(mat) = names(envVariables[[h]])
										# buffer[[t]] = mat
										mat
									}		
								for (t in 1:length(buffer))
									{
										buffer[[t]][!is.finite(buffer[[t]][])] = NA
										buffer[[t]][buffer[[t]][]==-1] = NA
										buffer[[t]][buffer[[t]][]==-777] = NA
										distancesSim[[t]][,h] = buffer[[t]][]
									}	
							}
					}
				if (impactOnDispersalLocation == TRUE)
					{
						simMeanEnvValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simMeanEnvDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simRateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						for (h in 2:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										envValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
										for (i in 1:length(ancestralNodes))
											{
												ancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
												envValues = envValues + raster::extract(envVariables[[h]], cbind(fromCoorRand[[t]][ancestralBranch,1],fromCoorRand[[t]][ancestralBranch,2]))
											}
										simMeanEnvValues[t,h] = mean(c(envValues,raster::extract(envVariables[[h]], cbind(toCoorRand[[t]][,1],toCoorRand[[t]][,2]))), na.rm=T)
										diffs1 = raster::extract(envVariables[[h]], fromCoorRand[[t]])-raster::extract(envVariables[[h]], toCoorRand[[t]])
										simRateOfPositiveDifferences[t,h] = sum(diffs1[!is.na(diffs1)] > 0)/length(diffs1[!is.na(diffs1)])
										diffs2 = raster::extract(envVariables[[h]], toCoorRand[[t]])-raster::extract(envVariables[[h]], fromCoorRand[[t]])
										simMeanEnvDifferences[t,h] = mean(diffs2, na.rm=T)
									}
							}
						colNames = c()
						for (h in 1:length(envVariables)) colNames = cbind(colNames, names(envVariables[[h]]))
						buffer = simMeanEnvValues; colnames(buffer) = colNames
						write.csv(buffer, paste(outputName,"_mean_evironmental_values_at_node_locations_randomisation_",s,".csv",sep=""), row.names=F, quote=F)
					}

				# 3.2. Linear regressions analyses on randomised datasets
				if (impactOnDiffusionVelocity == TRUE)
					{
						simUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						for (h in 1:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										durationTimes = dispersalTime[[t]]; distances1 = distancesSim[[t]][,h]
										LM1 = lm(as.formula("durationTimes ~ distances1"))
										simUniLRRsquares1[t,h] = summary(LM1)$r.squared
										simUniLRDeltaRsquares1[t,h] = simUniLRRsquares1[t,h]-simUniLRRsquares1[t,1]
										durationTimes4 = 4*dispersalTime[[t]]; squareDistances = distancesSim[[t]][,h]^2
										LM2 = lm(as.formula("durationTimes4 ~ squareDistances"))
										simUniLRRsquares2[t,h] = summary(LM2)$r.squared
										simUniLRDeltaRsquares2[t,h] = simUniLRRsquares2[t,h]-simUniLRRsquares2[t,1]
									}
							}
					}
				if (impactOnDiffusionVelocity == TRUE) H = 2
				if (impactOnDispersalLocation == TRUE) H = 2
					
				# 3.3. Computations of the BF values for this randomisation	
				if (impactOnDiffusionVelocity == TRUE)
					{
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
					}
				if (impactOnDispersalLocation == TRUE)
					{
						for (h in H:length(envVariables))
							{
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simMeanEnvValues[t,h]))&(!is.na(meanEnvValues[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (meanEnvValues[t,h] < simMeanEnvValues[t,h]) c = c+1
													}	else	{
														if (meanEnvValues[t,h] > simMeanEnvValues[t,h]) c = c+1
													}
											}	else	{
												missingValues = missingValues + 1
											}
									}
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								meanEnvValuesRandomisationBFs[h,s] = round(bf, 4)
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simMeanEnvDifferences[t,h]))&(!is.na(meanEnvDifferences[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (meanEnvDifferences[t,h] < simMeanEnvDifferences[t,h]) c = c+1
													}	else	{
														if (meanEnvDifferences[t,h] > simMeanEnvDifferences[t,h]) c = c+1
													}
											}	else	{
												missingValues = missingValues + 1
											}
									}
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								meanEnvDifferencesRandomisationBFs[h,s] = round(bf, 4)				
								c = 0; missingValues = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if ((!is.na(simRateOfPositiveDifferences[t,h]))&(!is.na(rateOfPositiveDifferences[t,h])))
											{
												if (resistances[h] == TRUE)
													{
														if (rateOfPositiveDifferences[t,h] > simRateOfPositiveDifferences[t,h]) c = c+1
													}	else	{
														if (rateOfPositiveDifferences[t,h] < simRateOfPositiveDifferences[t,h]) c = c+1
													}
											}	else	{
												missingValues = missingValues + 1
											}
									}	
								f = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
								rateOfPositiveDifferencesRandomisationBFs[h,s] = round(bf, 4)
							}
					}
			} # End of randomisations

		# 3.4. Creation of the text files containing the randomisation results
		if (nberOfExtractionFiles > 9)
			{
				envVariablesNames = c()
				for (h in 1:length(envVariables))
					{
						envVariablesNames = c(envVariablesNames, names(envVariables[[h]]))
					}
				if (impactOnDiffusionVelocity == TRUE)
					{
						colNames1 = c(); colNames2 = c()
						for (s in 1:nberOfRandomisations)
							{
								colNames1 = c(colNames1, paste("BFs_Q_LR1_randomisation_",s,sep=""))
								colNames2 = c(colNames2, paste("BFs_Q_LR2_randomisation_",s,sep=""))
							}
						row.names(uniLRdeltaRsquaresRandomisationBFs1) = envVariablesNames
						colnames(uniLRdeltaRsquaresRandomisationBFs1) = colNames1
						row.names(uniLRdeltaRsquaresRandomisationBFs2) = envVariablesNames
						colnames(uniLRdeltaRsquaresRandomisationBFs2) = colNames2
						if (reportingBothQstats == TRUE)
							{
								buffer = cbind(uniLRdeltaRsquaresRandomisationBFs1,uniLRdeltaRsquaresRandomisationBFs2)
							}	else	{
								buffer = uniLRdeltaRsquaresRandomisationBFs1
							}
						write.table(buffer, paste(outputName,"_randomisation_BF_results.txt",sep=""), quote=F, sep="\t")
					}
				if (impactOnDispersalLocation == TRUE)
					{
						row.names(meanEnvValuesRandomisationBFs) = envVariablesNames
						colnames(meanEnvValuesRandomisationBFs) = "BF"
						fileName = paste(outputName,"_position_E_BF_results.txt",sep="")
						write.table(meanEnvValuesRandomisationBFs, fileName, quote=F, sep="\t")
						row.names(meanEnvDifferencesRandomisationBFs) = envVariablesNames
						colnames(meanEnvDifferencesRandomisationBFs) = "BF"
						fileName = paste(outputName,"_direction_R_BF_results.txt",sep="")
						write.table(meanEnvDifferencesRandomisationBFs, fileName, quote=F, sep="\t")
					}
			}
	}
}
