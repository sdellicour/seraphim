isolationByResistance = function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), pathModel=2, resistances=list(), avgResistances=list(), 
								 fourCells=FALSE, nberOfRandomisations=0, randomProcedure=3, outputName="", showingPlots=FALSE, nberOfCores=1, 
								 OS="Unix", juliaCSImplementation=FALSE, simulations=FALSE, randomisations=FALSE, minimumConvexHull=TRUE) {

reportingBothQstats = TRUE; # simulations = FALSE; randomisations = FALSE
impactOnIsolation = TRUE
# registerDoMC(cores=nberOfCores)
nberOfCores_CS = 1
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
if (pathModel == 1) straightLineDistance = TRUE
if (pathModel == 2) leastCostDistance = TRUE
if (pathModel == 3) circuitscapeDistance = TRUE
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
nberOfTipNodes = rep(NA, nberOfExtractionFiles)
nberOfConnections = rep(NA, nberOfExtractionFiles); totalNberOfConnections = 0
node1 = list(); node2 = list(); fromCoor = list(); toCoor = list(); samplingCoor = list()
dispersalTime = list(); treeIDs = list()
if (impactOnIsolation == TRUE)
	{
		distances = list()
	}
for (t in 1:nberOfExtractionFiles)
	{
		if (nchar(localTreesDirectory) == 0)
			{
				fileName = paste(extractionFileName,"_",t,".csv",sep="")
			}	else	{
				fileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
			}	
		data = read.csv(fileName, head=T); nberOfConnections[t] = dim(data)[1]
		node1[[t]] = matrix(nrow=dim(data)[1], ncol=1); node1[[t]][] = data[,"node1"]
		node2[[t]] = matrix(nrow=dim(data)[1], ncol=1); node2[[t]][] = data[,"node2"]
		fromCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); fromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
		toCoor[[t]] = matrix(nrow=dim(data)[1], ncol=2); toCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
		tipIndices = which(!data[,"node2"]%in%data[,"node1"]); samplingCoor[[t]] = matrix(nrow=length(tipIndices), ncol=2)
		samplingCoor[[t]][] = cbind(data[tipIndices,"endLon"], data[tipIndices,"endLat"]); nberOfTipNodes[t] = length(tipIndices)
		dispersalTime[[t]] = matrix(nrow=dim(data)[1], ncol=1); colnames(dispersalTime[[t]]) = "dispersalTime"
		dispersalTime[[t]][] = data[,"endYear"]-data[,"startYear"]
		totalNberOfConnections = totalNberOfConnections + dim(data)[1]
		if (impactOnIsolation == TRUE)
			{
				nRows = ((nberOfTipNodes[t]*nberOfTipNodes[t])-nberOfTipNodes[t])/2
				distances[[t]] = matrix(nrow=nRows, ncol=length(envVariables))
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

# 1. Computation of all the patristic and environmental distances

cat("Computing patristic distances","\n",sep="")
patristicDistances = list(); buffer = list()
# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
for (t in 1:nberOfExtractionFiles) {
		if (nchar(localTreesDirectory) == 0)
			{
				fileName = paste(extractionFileName,"_",t,".csv",sep="")
			}	else	{
				fileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
			}	
		data = read.csv(fileName, head=T); nberOfConnections[t] = dim(data)[1]
		tipNodeIndices = which(!data[,"node2"]%in%data[,"node1"])
		distTree = matrix(nrow=length(tipNodeIndices), ncol=length(tipNodeIndices))
		for (i in 2:dim(distTree)[1])
			{
				for (j in 1:(i-1))
					{
						index1 = tipNodeIndices[i]
						index2 = tipNodeIndices[j]
						indices1 = index1; root = FALSE
						while (root == FALSE)
							{	
								if (data[indices1[length(indices1)],"node1"]%in%data[,"node2"])
									{
										indices1 = c(indices1, which(data[,"node2"]==data[indices1[length(indices1)],"node1"]))
									}	else	{
										root = TRUE
									}
							}
						indices2 = index2; root = FALSE
						while (root == FALSE)
							{	
								if (data[indices2[length(indices2)],"node1"]%in%data[,"node2"])
									{
										indices2 = c(indices2, which(data[,"node2"]==data[indices2[length(indices2)],"node1"]))
									}	else	{
										root = TRUE
									}
							}
						indices3 = indices1[which(indices1%in%indices2)]; patristic_dis = NULL
						if (length(indices3) == 0)
							{
								patristic_dis = sum(data[c(indices1,indices2),"length"])
							}	else	{
								patristic_dis = sum(data[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
							}
						distTree[i,j] = patristic_dis; distTree[i,j] = patristic_dis
					}
			}
		patristicDistances[[t]] = distTree[lower.tri(distTree)]
		# distTree[lower.tri(distTree)]
	}
for (t in 1:length(buffer))
	{
		patristicDistances[[t]] = buffer[[t]]
	}
for (h in 1:length(envVariables))
	{
		if (impactOnIsolation == TRUE)
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
				# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
				for (t in 1:nberOfExtractionFiles) {
						mat = matrix(nrow=dim(samplingCoor[[t]])[1], ncol=dim(samplingCoor[[t]])[1])
						if (straightLineDistance == TRUE)
							{
								linesList = list(); c = 0
								for (i in 2:dim(samplingCoor[[t]])[1])
									{
										for (j in 1:(i-1))
											{
												points = rbind(samplingCoor[[t]][i,], samplingCoor[[t]][j,])
												c = c+1; linesList[[c]] = Lines(list(Line(points)),c)
											}
									}
								lines = SpatialLines(linesList)
								extractions = raster::extract(hullRasters[[h]], lines); c = 0
								for (i in 2:length(samplingCoor[[t]][,1]))
									{
										for (j in 1:(i-1))
											{
												c = c+1; mat[i,j] = sum(extractions[[c]], na.rm=T); mat[j,i] = mat[i,j]
											}
									}
							}
						if (leastCostDistance == TRUE)
							{
								mat[] = costDistance(trEnvVariableCorr, samplingCoor[[t]], samplingCoor[[t]])
							}
						if (circuitscapeDistance == TRUE)
							{
								envVariableName = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
								samplingPointNotNA = which(!(is.na(raster::extract(hullRasters[[h]], samplingCoor[[t]][]))))
								samplingCoor_temp = samplingCoor[[t]][samplingPointNotNA,]
								if (juliaCSImplementation == FALSE)
									{
										mat[samplingPointNotNA,samplingPointNotNA] = circuitScape1(hullRasters[[h]], envVariableName, resistances[[h]], avgResistances[[h]],
																			 					   fourCells, samplingCoor_temp, samplingCoor_temp, OS, outputName, t, nberOfCores_CS)
									}	else	{
										mat[samplingPointNotNA,samplingPointNotNA] = circuitScape2(hullRasters[[h]], envVariableName, resistances[[h]], avgResistances[[h]],
																			 					   fourCells, samplingCoor_temp, samplingCoor_temp, OS, outputName, t, nberOfCores_CS)
									}
							}
						buffer[[t]] = mat
						# mat
					}		
				for (t in 1:length(buffer))
					{
						for (i in 1:dim(buffer[[t]])[1])
							{
								for (j in 1:dim(buffer[[t]])[2])
									{
										if (!is.finite(buffer[[t]][i,j])) buffer[[t]][i,j] = NA # least-cost case
										if ((!is.na(buffer[[t]][i,j]))&&(buffer[[t]][i,j] == -1)) buffer[[t]][i,j] = NA # CircuitScape case
										if ((!is.na(buffer[[t]][i,j]))&&(buffer[[t]][i,j] == -777)) buffer[[t]][i,j] = NA # NA value in Circuitscape
									}
							}
						distances[[t]][,h] = buffer[[t]][lower.tri(buffer[[t]])]
					}
			}
	}
envVariableNames = names(envVariables[[1]])
for (h in 2:length(envVariables))
	{
		envVariableNames = cbind(envVariableNames, names(envVariables[[h]]))
	}
if ((nberOfExtractionFiles == 1)&(impactOnIsolation == TRUE))
	{
		mat = cbind(patristicDistances[[1]][], distances[[1]][,1:length(envVariables)])
		columnNames = "patristic_distances"
		for (h in 1:length(envVariables))
			{
				columnNames = cbind(columnNames, names(envVariables[[h]]))
			}
		colnames(mat) = columnNames
		write.table(mat, file=paste(outputName,"_env_distances.txt",sep=""), row.names=F, quote=F, sep="\t")
	}
if ((file.exists(outputName))&(impactOnIsolation == TRUE))
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

if (impactOnIsolation == TRUE)
	{
		realUniLRcoefficients = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniLRRsquares = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realUniDeltaRsquares = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realRP2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
		realDeltaRP2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables)); colNames = c()
		for (t in 1:nberOfExtractionFiles) colNames = c(colNames, paste("tree_",treeIDs[[t]],sep=""))
		for (h in 1:length(envVariables))
			{
				nRowsMax = length(patristicDistances[[1]])
				if (nberOfExtractionFiles > 1)
					{
						for (t in 2:nberOfExtractionFiles)
							{
								if (nRowsMax < length(patristicDistances[[t]])) nRowsMax = length(patristicDistances[[t]])
							}
					}
				for (t in 1:nberOfExtractionFiles)
					{
						patristicDistance = patristicDistances[[t]]; environmentalDistance = distances[[t]][,h]
						patristicDistance = patristicDistance[which(!is.na(environmentalDistance))]
						environmentalDistance = environmentalDistance[which(!is.na(environmentalDistance))]
						environmentalDistanceLog = log(environmentalDistance+1)
						LM = lm(as.formula("patristicDistance ~ environmentalDistanceLog"))
						realUniLRcoefficients[t,h] = summary(LM)$coefficients[2,"Estimate"]
						realUniLRRsquares[t,h] = summary(LM)$r.squared
						realRP2[t,h] = cor(patristicDistance,log(environmentalDistance+1), method="pearson")
					}
			}
		for (h in 2:length(envVariables))
			{
				realUniDeltaRsquares[,h] = (realUniLRRsquares[,h]-realUniLRRsquares[,1])
				realDeltaRP2[,h] = realRP2[,h]-realRP2[,1]
			}
		if (nberOfExtractionFiles > 1)
			{
				uniLRcoefficientsNames = c(); uniLRRsquaresNames = c(); uniLRdeltaRsquaresNames = c(); realRP2Names = c(); realDeltaRP2Names = c()
				mat1 = cbind(realUniLRcoefficients[,1:dim(realUniLRcoefficients)[2]], realUniLRRsquares[,1:dim(realUniLRRsquares)[2]], realUniDeltaRsquares[,2:dim(realUniDeltaRsquares)[2]])
				mat2 = cbind(realRP2[,1:dim(realRP2)[2]], realDeltaRP2[,2:dim(realDeltaRP2)[2]])
				for (h in 1:length(envVariables))
					{
						uniLRcoefficientsNames = cbind(uniLRcoefficientsNames, paste("LR_coefficients_",names(envVariables[[h]]),sep=""))	
						uniLRRsquaresNames = cbind(uniLRRsquaresNames, paste("LR_R2_",names(envVariables[[h]]), sep=""))				
						if (h > 1) uniLRdeltaRsquaresNames = cbind(uniLRdeltaRsquaresNames,paste("LR_Q_",names(envVariables[[h]]),sep=""))
						realRP2Names = cbind(realRP2Names, paste("Pearson_correlation_",names(envVariables[[h]]),sep=""))	
						if (h > 1) realDeltaRP2Names = cbind(realDeltaRP2Names, paste("Pearson_correlation_difference_",names(envVariables[[h]]),sep=""))
					}
				names1 = cbind(uniLRcoefficientsNames, uniLRRsquaresNames, uniLRdeltaRsquaresNames); colnames(mat1) = names1
				names2 = cbind(realRP2Names, realDeltaRP2Names); colnames(mat2) = names2
				write.table(mat1, file=paste(outputName,"_linear_regression_results.txt",sep=""), row.names=F, quote=F, sep="\t")
				# write.table(mat2, file=paste(outputName,"_Pearson_correlation_results.txt",sep=""), row.names=F, quote=F, sep="\t")
			}
	}

# 3. Randomisation step (randomisation of phylogenetic branch positions)

if (nberOfRandomisations > 0)
	{	
		if (impactOnIsolation == TRUE)
			{
				uniLRdeltaRsquaresRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
				deltaRP2RandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
			}
		for (s in 1:nberOfRandomisations)
			{
				if (impactOnIsolation == TRUE)
					{
						uniLRdeltaRsquaresRand = matrix(0, nrow=nberOfExtractionFiles, ncol=nberOfRandomisations)
						distancesSim = list()
						for (t in 1:nberOfExtractionFiles)
							{
								nRows = ((nberOfTipNodes[t]*nberOfTipNodes[t])-nberOfTipNodes[t])/2
								distancesSim[[t]] = matrix(nrow=nRows, ncol=length(envVariables))
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
						samplingCoor = list(); dispersalTime = list(); treeIDs = list()
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
								tipIndices = which(!data[,"node2"]%in%data[,"node1"]); samplingCoor[[t]] = matrix(nrow=length(tipIndices), ncol=2)
								samplingCoor[[t]][] = cbind(data[tipIndices,"endLon"], data[tipIndices,"endLat"])
								dispersalTime[[t]] = matrix(nrow=dim(data)[1], ncol=1); colnames(dispersalTime[[t]]) = "dispersalTime"
								dispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
								if (("treeID"%in%colnames(data)) == TRUE)
									{
										treeIDs[[t]] = data[1,"treeID"]
									}	else	{
										treeIDs[[t]] = "noTreeID"
									}	
							}
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor
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
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor				
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
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor
					}
				if (branchRandomisation3 == TRUE)
					{
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor
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
								samplingCoorRand[[t]] = toCoorRand[[t]][which(!node2[[t]]%in%node1[[t]]),]
							}
					}
				if (branchRandomisation2 == TRUE)
					{
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor
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
								samplingCoorRand[[t]] = toCoorRand[[t]][which(!node2[[t]]%in%node1[[t]]),]
							}
					}			
				if (branchRandomisation1 == TRUE)
					{
						simRasters = list(); simRasters = hullRasters
						cat("Analysis of randomised branch positions ",s,"\n",sep="")
						fromCoorRand = fromCoor; toCoorRand = toCoor; samplingCoorRand = samplingCoor
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
								samplingCoorRand[[t]] = toCoorRand[[t]][which(!node2[[t]]%in%node1[[t]]),]
							}
					}		
				if (impactOnIsolation == TRUE)
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
								# buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
								for (t in 1:nberOfExtractionFiles) {
										mat = matrix(nrow=dim(samplingCoorRand[[t]])[1], ncol=dim(samplingCoorRand[[t]])[1])
										if (straightLineDistance == TRUE)
											{	
												linesList = list(); c = 0
												for (i in 2:dim(samplingCoorRand[[t]])[1])
													{
														for (j in 1:(i-1))
															{
																points = rbind(samplingCoorRand[[t]][i,], samplingCoorRand[[t]][j,])
																c = c+1; linesList[[c]] = Lines(list(Line(points)),c)
															}
													}
												lines = SpatialLines(linesList)
												extractions = raster::extract(hullRasters[[h]], lines); c = 0
												for (i in 2:length(samplingCoorRand[[t]][,1]))
													{
														for (j in 1:(i-1))
															{
																c = c+1; mat[i,j] = sum(extractions[[c]], na.rm=T); mat[j,i] = mat[i,j]
															}
													}
											}	
										if (leastCostDistance == TRUE)
											{	
												mat[] = costDistance(trEnvVariableCorr, samplingCoorRand[[t]], samplingCoorRand[[t]])
											}
										if (circuitscapeDistance == TRUE)
											{
												simRasterName = paste("CS_rasters/",names(simRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
												samplingPointNotNA = which(!(is.na(raster::extract(simRasters[[h]], samplingCoorRand[[t]][]))))
												samplingCoor_temp = samplingCoorRand[[t]][samplingPointNotNA,]
												if (juliaCSImplementation == FALSE)
													{
														mat[samplingPointNotNA,samplingPointNotNA] = circuitScape1(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], 
																												   fourCells, samplingCoor_temp, samplingCoor_temp, OS, outputName, t, 
																												   nberOfCores_CS)
													}	else	{
														mat[samplingPointNotNA,samplingPointNotNA] = circuitScape2(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], 
																												   fourCells, samplingCoor_temp, samplingCoor_temp, OS, outputName, t, 
																												   nberOfCores_CS)
													}
											}
										buffer[[t]] = mat
										# mat
									}		
								for (t in 1:length(buffer))
									{
										for (i in 1:dim(buffer[[t]])[1])
											{
												for (j in 1:dim(buffer[[t]])[2])
													{
														if (!is.finite(buffer[[t]][i,j])) buffer[[t]][i,j] = NA # least-cost case
														if ((!is.na(buffer[[t]][i,j]))&&(buffer[[t]][i,j] == -1)) buffer[[t]][i,j] = NA # CircuitScape case
														if ((!is.na(buffer[[t]][i,j]))&&(buffer[[t]][i,j] == -777)) buffer[[t]][i,j] = NA # NA value in Circuitscape
													}
											}
										distancesSim[[t]][,h] = buffer[[t]][lower.tri(buffer[[t]])]
									}
							}
					}

				# 3.2. Linear regressions analyses on randomised datasets
				if (impactOnIsolation == TRUE)
					{
						simUniLRRsquares = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simUniLRDeltaRsquares = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simRP2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						simDeltaRP2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
						for (h in 1:length(envVariables))
							{
								for (t in 1:nberOfExtractionFiles)
									{
										patristicDistance = patristicDistances[[t]]; environmentalDistance = distancesSim[[t]][,h]
										patristicDistance = patristicDistance[which(!is.na(environmentalDistance))]
										environmentalDistance = environmentalDistance[which(!is.na(environmentalDistance))]
										environmentalDistanceLog = log(environmentalDistance+1)
										LM = lm(as.formula("patristicDistance ~ environmentalDistanceLog"))
										simUniLRRsquares[t,h] = summary(LM)$r.squared
										simUniLRDeltaRsquares[t,h] = simUniLRRsquares[t,h]-simUniLRRsquares[t,1]
										simRP2[t,h] = cor(patristicDistance,log(environmentalDistance+1), method="pearson")
										simDeltaRP2[t,h] = simRP2[t,h]-simRP2[t,1]
									}
							}
					}
				if (impactOnIsolation == TRUE) H = 2
					
				# 3.3. Computations of the BF values for this randomisation	
				if (impactOnIsolation == TRUE)
					{
						for (h in H:length(envVariables))
							{
								c = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if (simUniLRDeltaRsquares[t,h] < realUniDeltaRsquares[t,h]) c = c+1
									}
								f = c/nberOfExtractionFiles; bf = f/(1-f)
								uniLRdeltaRsquaresRandomisationBFs[h,s] = round(bf, 4)
								c = 0
								for (t in 1:nberOfExtractionFiles)
									{
										if (simDeltaRP2[t,h] < realDeltaRP2[t,h]) c = c+1
									}
								f = c/nberOfExtractionFiles; bf = f/(1-f)
								# deltaRP2RandomisationBFs[h,s] = round(bf, 4)
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
				if (impactOnIsolation == TRUE)
					{
						colNames = c()
						for (s in 1:nberOfRandomisations)
							{
								colNames = c(colNames, paste("BFs_Q_LR_randomisation_",s,sep=""))
								# colNames = c(colNames, paste("BFs_rP2_randomisation_",s,sep=""))
							}
						row.names(uniLRdeltaRsquaresRandomisationBFs) = envVariablesNames
						colnames(uniLRdeltaRsquaresRandomisationBFs) = colNames
						write.table(uniLRdeltaRsquaresRandomisationBFs, paste(outputName,"_randomisation_BF_results.txt",sep=""), quote=F, sep="\t")
					}
			}
	}
}
