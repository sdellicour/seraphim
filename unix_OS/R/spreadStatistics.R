spreadStatistics <-
function(localTreesDirectory="", nberOfExtractionFiles=1, timeSlices=200, onlyTipBranches=F, showingPlots=TRUE, outputName=gsub(" ","_",date()), nberOfCores=1, slidingWindow=NA, simulations=FALSE, discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot=FALSE) {

\tnberOfStatistics = 6; treeIDs = c()
\tregisterDoMC(cores=nberOfCores)
  \t\t# (1) mean branch dispersal velocity
  \t\t# (2) weighted branch dispersal velocity
  \t\t# (3) branch dispersal velocity variation among branches (CV)
  \t\t# (4) original diffuson coefficient (Pybus et al. 2012)
  \t\t# (5) weighted diffuson coefficient (Trovao et al. 2015)
  \t\t# (6) diffusion coefficient variation among branches (CV)
\tmeanStatistics = matrix(nrow=(nberOfExtractionFiles), ncol=nberOfStatistics)
\tbranchVelocities = c() # not used, just to obtain an overall distribution of velocities
\tsd_var_velocity = matrix(nrow=(nberOfExtractionFiles), ncol=2) # not used either
\tmedianMeanStatistics = matrix(nrow=1, ncol=nberOfStatistics)
\tciMeanStatistics = matrix(nrow=2, ncol=nberOfStatistics)
\twaveFrontDistances1List = list()
\twaveFrontDistances2List = list()
\tmeanDispersalVelocityList = list()
\tweightedDispersalVelocityList = list()
\tnumberOfBranchesList = list()
\tdispersalOrientationList = list()
\tcat("Estimation of summary statistics", "\\n", sep="")
\tonlyOneAncestor = TRUE; extractionsWithMoreThanOneAncestors = c()
\tif ((timeSlices==0)|is.na(timeSlices)|is.null(timeSlices)) onlyOneAncestor = FALSE
\tdispersalVelocityGraph = FALSE
\tif (!is.na(slidingWindow)) dispersalVelocityGraph = TRUE
\tfor (t in 1:nberOfExtractionFiles)
\t\t{
\t\t\tif (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
\t\t\tif (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
\t\t\tdata = read.csv(fileName, h=T)
\t\t\tdata = data[with(data, order(endYear,startYear)),]
\t\t\tif (sum(!data[,"node1"]%in%data[,"node2"]) > 2)
\t\t\t\t{
\t\t\t\t\tonlyOneAncestor = FALSE; extractionsWithMoreThanOneAncestors = c(extractionsWithMoreThanOneAncestors, t)
\t\t\t\t}
\t\t\ttreeIDs = c(treeIDs, data[1,"treeID"])
\t\t\tif (onlyTipBranches == TRUE)
\t\t\t\t{
\t\t\t\t\tindices = c()
\t\t\t\t\tfor (i in 1:dim(data)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tn = data[i,"node2"]
\t\t\t\t\t\t\tif (length(data[data[,"node1"]==n,"node1"]) == 0)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tindices = c(indices, i)
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tdata = data[indices,]\t
\t\t\t\t}
\t\t\tnberOfConnections = dim(data)[1]
\t\t\tif (t == 1)
\t\t\t\t{
\t\t\t\t\tminLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
\t\t\t\t\tmaxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
\t\t\t\t\tminLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
\t\t\t\t\tmaxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
\t\t\t\t\tminStartYear = min(data[,"startYear"])
\t\t    \t\tminEndYear = min(data[,"startYear"])
\t\t    \t\tmaxEndYear = max(data[,"endYear"])
\t\t\t\t}\telse\t{
\t\t\t\t\tif (minLon > min(min(data[,"endLon"]),min(data[,"startLon"]))) minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
\t\t\t\t\tif (maxLon < max(max(data[,"endLon"]),max(data[,"startLon"]))) maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
\t\t\t\t\tif (minLat > min(min(data[,"endLat"]),min(data[,"startLat"]))) minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
\t\t\t\t\tif (maxLat < max(max(data[,"endLat"]),max(data[,"startLat"]))) maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
\t\t\t\t\tif (minStartYear > min(data[,"startYear"])) minStartYear = min(data[,"startYear"])
\t\t\t\t\tif (minEndYear > min(data[,"endYear"])) minEndYear = min(data[,"endYear"])
\t\t\t\t\tif (maxEndYear < max(data[,"endYear"])) maxEndYear = max(data[,"endYear"])
\t\t\t\t}
\t\t}
\tif ((onlyOneAncestor == FALSE)&(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE))
\t\t{
\t\t\tcat("Discarding",length(extractionsWithMoreThanOneAncestors),"extraction tables with more without a single","\\n")
\t\t\tcat("\\t","common ancestor for generating wavefront plots","\\n")
\t\t}
\twavefrontDistanceSlices = timeSlices
\twavefrontDistanceTimeInterval = (maxEndYear-minStartYear)/wavefrontDistanceSlices
\tdispersalVelocitySlices = timeSlices
\tdispersalVelocityTimeInterval = (maxEndYear-minStartYear)/dispersalVelocitySlices
\txLim = c(minStartYear, maxEndYear)
\tfor (t in 1:nberOfExtractionFiles)
\t\t{
\t\t\tif (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
\t\t\tif (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
\t\t\tdata = read.csv(fileName, h=T)
\t\t\tdata = data[with(data, order(endYear,startYear)),]
\t\t\tnberOfConnections = dim(data)[1]
\t\t\tdistances = matrix(0, nrow=dim(data)[1], ncol=1)
\t \t\tcolnames(distances) = c("greatCircleDist_km")
\t\t\tfor (i in 1:nberOfConnections)
\t\t\t\t{
\t\t\t\t\tx1 = cbind(data[i,"startLon"],data[i,"startLat"])
\t \t\t\t\tx2 = cbind(data[i,"endLon"],data[i,"endLat"])
\t \t\t\t\tdistances[i] = rdist.earth(x1, x2, miles=F, R=NULL)
\t\t\t\t}
\t\t\tif (length(grep("greatCircleDist_km",colnames(data))) > 0)
\t\t\t\t{
\t\t\t\t\tdata[,"greatCircleDist_km"] = distances
\t\t\t\t}\telse\t{
\t\t\t\t\tdata = cbind(data,distances)
\t\t\t\t}
\t\t\tbranchMeasures = matrix(nrow=nberOfConnections, ncol=2)
\t\t\tweightedDispersalVelocity_numerator = 0; weightedDispersalVelocity_denominator = 0
\t\t\tweightedDiffusionCoefficient_numerator = 0; weightedDiffusionCoefficient_denominator = 0
\t\t\tfor (i in 1:nberOfConnections)
\t\t\t\t{
\t\t\t\t\tdispersalTime = data[i,"endYear"]-data[i,"startYear"]
\t\t    \t\t\tbranchVelocity = data[i,"greatCircleDist_km"]/(dispersalTime)
\t\t    \t\t\tbranchOriginalDiffusionCoefficient = (data[i,"greatCircleDist_km"]^2)/(4*dispersalTime)
\t\t    \t\t\tbranchMeasures[i,1] = branchVelocity
\t\t    \t\t\tbranchMeasures[i,2] = branchOriginalDiffusionCoefficient
\t\t    \t\t\tweightedDispersalVelocity_numerator = weightedDispersalVelocity_numerator + data[i,"greatCircleDist_km"]
\t\t    \t\t\tweightedDispersalVelocity_denominator = weightedDispersalVelocity_denominator + dispersalTime
\t\t    \t\t\tweightedDiffusionCoefficient_numerator = weightedDiffusionCoefficient_numerator + (data[i,"greatCircleDist_km"]^2)
\t\t    \t\t\tweightedDiffusionCoefficient_denominator = weightedDiffusionCoefficient_denominator + (4*dispersalTime)
\t\t  \t\t}
\t\t  \tbranchVelocities = c(branchVelocities, branchMeasures[,1])
\t\t\tsd_var_velocity[t,] = cbind(sd(branchMeasures[,1]), var(branchMeasures[,1]))
\t\t\tmeanStatistics[t,1] = mean(branchMeasures[,1])
\t\t\tmeanStatistics[t,2] = weightedDispersalVelocity_numerator/weightedDispersalVelocity_denominator
\t\t\tmeanStatistics[t,3] = sd(branchMeasures[,1])/mean(branchMeasures[,1])
\t\t\tmeanStatistics[t,4] = mean(branchMeasures[,2])
\t\t    meanStatistics[t,5] = weightedDiffusionCoefficient_numerator/weightedDiffusionCoefficient_denominator
\t\t\tmeanStatistics[t,6] = sd(branchMeasures[,2])/mean(branchMeasures[,2])
\t\t}
\tif ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&((onlyOneAncestor == TRUE)|(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE)))
\t\t{
\t\t\tbuffer = list()
\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {\t\t\t\t
\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
\t\t\t\t\t\t\tif (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
\t\t\t\t\t\t\tdata = read.csv(fileName, h=T)
\t\t\t\t\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\t\t\t\t\tnberOfConnections = dim(data)[1]
\t\t\t\t\t\t\twaveFrontDistances1 = matrix(nrow=wavefrontDistanceSlices+1, ncol=2) # spatial wavefront distance
\t\t\t\t\t\t\twaveFrontDistances2 = matrix(nrow=wavefrontDistanceSlices+1, ncol=2) # patristic wavefront distance
\t\t\t\t\t\t\twaveFrontDistances1[,1] = seq(minStartYear,maxEndYear,wavefrontDistanceTimeInterval)
\t\t\t\t\t\t\twaveFrontDistances2[,1] = seq(minStartYear,maxEndYear,wavefrontDistanceTimeInterval)
\t\t\t\t\t\t\tdata = data[order(data[,"endYear"]),]
\t\t\t\t\t\t\tancestralIndex = which(data[,"startYear"]==min(data[,"startYear"]))[1]
\t\t\t\t\t\t\tancestralNode = data[ancestralIndex,"node1"]
\t\t\t\t\t\t\tstartingYear = data[ancestralIndex,"startYear"]
\t\t\t\t\t\t\toriginLocationX = data[ancestralIndex,"startLon"]
\t\t\t\t\t\t\toriginLocationY = data[ancestralIndex,"startLat"]
\t\t\t\t\t\t\toriginLocation = cbind(originLocationX, originLocationY)
\t\t\t\t\t\t\tmaxDistance1 = 0; maxDistance2 = 0
\t\t\t\t\t\t\tfor (i in 0:wavefrontDistanceSlices)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\ttime = minStartYear+(i*wavefrontDistanceTimeInterval)
\t\t\t\t\t\t\t\t\tfor (j in 1:nberOfConnections)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif ((data[j,"startYear"]<time)&(time<data[j,"endYear"]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\ttimeProportion = (time-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
\t\t\t\t\t\t\t\t\t\t\t\t\tpointLocationX = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion)
\t\t\t\t\t\t\t\t\t\t\t\t\tpointLocationY = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion)
\t\t\t\t\t\t\t\t\t\t\t\t\tpointLocation = cbind(pointLocationX, pointLocationY)
\t\t\t\t\t\t\t\t\t\t\t\t\tdistance1 = rdist.earth(originLocation, pointLocation, miles=F, R=NULL)\t\t\t\t    \t\t \t\t
\t\t\t\t\t\t    \t\t\t\t\t\t\tcoordinatesOfStartNode = cbind(data[j,"startLon"], data[j,"startLat"])
\t\t\t\t\t\t    \t\t\t\t\t\t\tdistance2 = rdist.earth(coordinatesOfStartNode, pointLocation, miles=F, R=NULL)
\t\t\t\t\t\t    \t\t\t\t\t\t\toriginNode = data[j,"node1"]
\t\t\t\t\t\t    \t\t\t\t\t\t\twhile (originNode != ancestralNode)
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t{
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t\tindex = which(data[,"node2"]==originNode)[1]
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t\tcoordinatesOfEndNode = cbind(data[index,"endLon"], data[index,"endLat"])
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t\tcoordinatesOfStartNode = cbind(data[index,"startLon"], data[index,"startLat"])
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t\tdistance2 = distance2 + rdist.earth(coordinatesOfStartNode, coordinatesOfEndNode, miles=F, R=NULL)
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t\toriginNode = data[index,"node1"]
\t\t\t\t\t\t    \t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tif (distance1 > maxDistance1) maxDistance1 = distance1
\t\t\t\t\t\t    \t\t\t\t\t\t\tif (distance2 > maxDistance2) maxDistance2 = distance2
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\twaveFrontDistances1[i+1,2] = maxDistance1
\t\t\t\t\t\t    \t\twaveFrontDistances2[i+1,2] = maxDistance2
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tcolnames(waveFrontDistances1) = c("endYear", "waveFrontDistances1")
\t\t\t\t\t\t\tcolnames(waveFrontDistances2) = c("endYear", "waveFrontDistances2")
\t\t\t\t\t\t\twaveFrontDistances = list()
\t\t\t\t\t\t\twaveFrontDistances[[1]] = waveFrontDistances1
\t\t\t\t\t\t\twaveFrontDistances[[2]] = waveFrontDistances2
\t\t\t\t\t\t\t# buffer[[t]] = waveFrontDistances
\t\t\t\t\t\t\twaveFrontDistances
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tfor (t in 1:length(buffer))
\t\t\t\t{
\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors)
\t\t\t\t\t\t{
\t\t\t\t\t\t\twaveFrontDistances1List[[t]] = buffer[[t]][[1]]
\t\t\t\t\t\t\twaveFrontDistances2List[[t]] = buffer[[t]][[2]]
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tmaxDistance1 = 0; maxDistance2 = 0
\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t{
\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (max(waveFrontDistances1List[[t]][,2]) > maxDistance1) maxDistance1 = max(waveFrontDistances1List[[t]][,2])
\t\t\t\t\t\t\tif (max(waveFrontDistances2List[[t]][,2]) > maxDistance2) maxDistance2 = max(waveFrontDistances2List[[t]][,2])
\t\t\t\t\t\t}
\t\t\t\t}
\t\t}
\tif ((onlyTipBranches == FALSE)&(dispersalVelocityGraph == TRUE)&(nberOfExtractionFiles > 1))
\t\t{
\t\t\tstartEndTimes = matrix(nrow=dispersalVelocitySlices, ncol=3)
\t\t\tfor (i in 1:dispersalVelocitySlices)
\t\t\t\t{
\t\t\t\t\ttime = minStartYear+((i-1)*dispersalVelocityTimeInterval)+(dispersalVelocityTimeInterval/2)
\t\t\t\t\tstartTime = time - (slidingWindow/2)
\t\t\t\t\tendTime = time + (slidingWindow/2)
\t\t\t\t\tstartEndTimes[i,1:3] = cbind(time, startTime, endTime)
\t\t\t\t}
\t\t\tbuffer = list()
\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {\t\t
\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\tif (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
\t\t\t\t\tif (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
\t\t\t\t\tdata = read.csv(fileName, h=T)
\t\t\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\t\t\tnberOfConnections = dim(data)[1]
\t\t\t\t\tmeanDispersalVelocities = matrix(nrow=dispersalVelocitySlices, ncol=2)
\t\t\t\t\tweightedDispersalVelocities = matrix(nrow=dispersalVelocitySlices, ncol=2)
\t\t\t\t\tnumberOfBranches = matrix(nrow=dispersalVelocitySlices, ncol=2)
\t\t\t\t\tdata = data[order(data[,"endYear"]),]
\t\t\t\t\tfor (i in 1:dispersalVelocitySlices)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tn = 0; vS = 0; dS = 0; tS = 0 
\t\t\t\t\t\t\ttime = startEndTimes[i,1]
\t\t\t\t\t\t\tstartTime = startEndTimes[i,2]
\t\t\t\t\t\t\tendTime = startEndTimes[i,3]
\t\t\t\t\t\t\tfor (j in 1:nberOfConnections)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tbranchInwavefrontDistanceTimeInterval = FALSE
\t\t\t\t\t\t\t\t\tif ((data[j,"startYear"]<startTime)&(data[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
\t\t\t\t\t\t\t\t\tif ((data[j,"startYear"]>startTime)&(data[j,"startYear"]<endTime)) branchInwavefrontDistanceTimeInterval = TRUE
\t\t\t\t\t\t\t\t\tif ((data[j,"endYear"]>startTime)&(data[j,"endYear"]<endTime)) branchInwavefrontDistanceTimeInterval = TRUE
\t\t\t\t\t\t\t\t\tif (branchInwavefrontDistanceTimeInterval == TRUE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (data[j,"startYear"] > startTime) time1 = data[j,"startYear"]
\t\t\t\t\t\t\t\t\t\t\tif (data[j,"startYear"] <= startTime) time1 = startTime
\t\t\t\t\t\t\t\t\t\t\tif (data[j,"endYear"] < endTime) time2 = data[j,"endYear"]
\t\t\t\t\t\t\t\t\t\t\tif (data[j,"endYear"] >= endTime) time2 = endTime
\t\t\t\t\t\t\t\t\t\t\tbranchTimeInInterval = time2 - time1
\t\t\t\t\t\t\t\t\t\t\ttimeProportion1 = (time1-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
\t\t\t\t\t\t\t\t\t\t\tpointLocation1X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion1)
\t\t\t\t\t\t\t\t\t\t\tpointLocation1Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion1)
\t\t\t\t\t\t\t\t\t\t\tpointLocation1 = cbind(pointLocation1X, pointLocation1Y)
\t\t\t\t\t\t\t\t\t\t\ttimeProportion2 = (time2-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
\t\t\t\t\t\t\t\t\t\t\tpointLocation2X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion2)
\t\t\t\t\t\t\t\t\t\t\tpointLocation2Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion2)
\t\t\t\t\t\t\t\t\t\t\tpointLocation2 = cbind(pointLocation2X, pointLocation2Y)
\t\t\t\t\t\t\t\t\t\t\tbranchDistInInterval = rdist.earth(pointLocation1, pointLocation2, miles=F, R=NULL)\t\t\t\t    \t\t \t\t
\t\t\t\t    \t\t\t\t\t\t\tn = n+1; vS = vS + (branchDistInInterval/branchTimeInInterval)
\t\t\t\t    \t\t\t\t\t\t\tdS = dS + branchDistInInterval; tS = tS + branchTimeInInterval
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tmeanDispersalVelocities[i,1] = time
\t\t\t\t\t\t\tmeanDispersalVelocities[i,2] = vS/n
\t\t\t\t\t\t\tweightedDispersalVelocities[i,1] = time
\t\t\t\t\t\t\tweightedDispersalVelocities[i,2] = dS/tS
\t\t\t\t\t\t\tnumberOfBranches[i,1] = time
\t\t\t\t\t\t\tnumberOfBranches[i,2] = n
\t\t\t\t\t\t}
\t\t\t\t\tcolnames(meanDispersalVelocities) = c("year","meanDispersalVelocity")
\t\t\t\t\tcolnames(weightedDispersalVelocities) = c("year","weightedDispersalVelocity")
\t\t\t\t\tdispersalVelocities = list()
\t\t\t\t\tdispersalVelocities[[1]] = meanDispersalVelocities
\t\t\t\t\tdispersalVelocities[[2]] = weightedDispersalVelocities
\t\t\t\t\tdispersalVelocities[[3]] = numberOfBranches
\t\t\t\t\t# buffer[[t]] = dispersalVelocities
\t\t\t\t\tdispersalVelocities
\t\t\t\t}
\t\t\tfor (t in 1:length(buffer))
\t\t\t\t{
\t\t\t\t\tmeanDispersalVelocityList[[t]] = buffer[[t]][[1]]
\t\t\t\t\tweightedDispersalVelocityList[[t]] = buffer[[t]][[2]]
\t\t\t\t\tnumberOfBranchesList[[t]] = buffer[[t]][[3]]
\t\t\t\t}
\t\t}
\tfor (i in 1:nberOfStatistics)
\t\t{
\t\t\tmedianMeanStatistics[1,i] = median(meanStatistics[,i], na.rm=T)
\t\t\tquantiles = quantile(meanStatistics[,i], probs=c(0.025,0.975), na.rm=T)
\t\t\tciMeanStatistics[1,i] = as.numeric(quantiles[1])
\t\t\tciMeanStatistics[2,i] = as.numeric(quantiles[2])
\t\t}
\tcat("Median value of mean branch dispersal velocity = ",medianMeanStatistics[1,1],"\\n\t95% HPD = [",ciMeanStatistics[1,1],", ",ciMeanStatistics[2,1],"]","\\n",sep="")
\tcat("Median value of weighted branch dispersal velocity = ",medianMeanStatistics[1,2],"\\n\t95% HPD = [",ciMeanStatistics[1,2],", ",ciMeanStatistics[2,2],"]","\\n",sep="")\t
\tcat("Median value of original diffusion coefficient = ",medianMeanStatistics[1,4],"\\n\t95% HPD = [",ciMeanStatistics[1,4],", ",ciMeanStatistics[2,4],"]","\\n",sep="")\t
\tcat("Median value of weighted diffusion coefficient = ",medianMeanStatistics[1,5],"\\n\t95% HPD = [",ciMeanStatistics[1,5],", ",ciMeanStatistics[2,5],"]","\\n",sep="")\t
\tcolnames(meanStatistics) = c("mean_branch_dispersal_velocity", "weighted_branch_dispersal_velocity", "branch_dispersal_velocity_variation_among_branches", 
\t"original_diffusion_coefficient", "weighted_diffusion_coefficient", "diffusion_coefficient_variation_among_branches")
\twrite.table(meanStatistics, file=paste(outputName,"_estimated_dispersal_statistics.txt",sep=""), quote=F, row.names=F, sep="\\t")

\tLWD = 0.2
\tif (nberOfExtractionFiles > 1)
\t\t{
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName, "_mean_branch_dispersal_velocity_variation.pdf",sep=""), width=5, height=5)
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
\t\t\t# H = Hpi(cbind(log10(meanStatistics[,1]), meanStatistics[,3]))
\t\t\t# kde = kde(cbind(log10(meanStatistics[,1]), meanStatistics[,3]), H=H)
\t\t\tH = Hpi(cbind(meanStatistics[,1], meanStatistics[,3]))
\t\t\tkde = kde(cbind(meanStatistics[,1], meanStatistics[,3]), H=H)
\t\t\ttext = "Kernel density estimates of mean branch dispersal velocity parameters"
\t\t\tcolours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
\t\t\txLab = "mean branch velocity"; yLab = "mean branch velocity variation among branches"
\t\t\tplot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30")
\t\t\tlegend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_variation.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()

\t\t\tif (showingPlots) dev.new(width=5,height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""),width=5,height=5)
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
\t\t\t# H = Hpi(cbind(log10(meanStatistics[,2]), meanStatistics[,3]))
\t\t\t# kde = kde(cbind(log10(meanStatistics[,2]), meanStatistics[,3]), H=H)
\t\t\tH = Hpi(cbind(meanStatistics[,2], meanStatistics[,3]))
\t\t\tkde = kde(cbind(meanStatistics[,2], meanStatistics[,3]), H=H)
\t\t\ttext = "Kernel density estimates of weighted branch dispersal velocity parameters"
\t\t\tcolours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
\t\t\txLab = "weighted dispersal velocity"; yLab="weighted dispersal velocity variation among branches"
\t\t\tplot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30")
\t\t\tlegend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()

\t\t\tif (showingPlots) dev.new(width=5,height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""),width=5,height=5)
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
\t\t\t# H = Hpi(cbind(log10(meanStatistics[,4]), meanStatistics[,6]))
\t\t\t# kde = kde(cbind(log10(meanStatistics[,4]), meanStatistics[,6]), H=H)
\t\t\tH = Hpi(cbind(meanStatistics[,4], meanStatistics[,6]))
\t\t\tkde = kde(cbind(meanStatistics[,4], meanStatistics[,6]), H=H)
\t\t\ttext1 = "Kernel density estimates of original diffusion coefficient parameters"; text2 = "(Pybus et al. 2012)"
\t\t\tcolours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
\t\t\txLab = "mean original diffusion coefficient"; yLab="diffusion coefficient variation among branches"
\t\t\tplot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text1, cex.main=0.6, col.main="gray30")
\t\t\tlegend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5,height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""), width=5, height=5)
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
\t\t\t# H = Hpi(cbind(log10(meanStatistics[,5]), meanStatistics[,6]))
\t\t\t# kde = kde(cbind(log10(meanStatistics[,5]), meanStatistics[,6]), H=H)
\t\t\tH = Hpi(cbind(meanStatistics[,5],meanStatistics[,6]))
\t\t\tkde = kde(cbind(meanStatistics[,5],meanStatistics[,6]),H=H)
\t\t\ttext1 = "Kernel density estimates of weighted diffusion coefficient parameters"; text2 = "(Trovao et al. 2015)"
\t\t\tcolours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
\t\t\txLab = "mean weighted diffusion coefficient"; yLab="diffusion coefficient variation among branches"
\t\t\tplot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text1, cex.main=0.6, col.main="gray30")
\t\t\tlegend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t}

\tif ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&((onlyOneAncestor == TRUE)|(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE)))
\t\t{
\t\t\tcat("Building wavefront distance evolution graphs", "\\n", sep="")
\t\t\twavefrontDistanceTimeInterval = (maxEndYear-minStartYear)/wavefrontDistanceSlices
\t\t\tslicedTimes = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\tlower_l_1 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\tupper_l_1 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\tlower_l_2 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\tupper_l_2 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances1Values = matrix(nrow=nberOfExtractionFiles-length(extractionsWithMoreThanOneAncestors), ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances1MeanValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances1MedianValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances2Values = matrix(nrow=nberOfExtractionFiles-length(extractionsWithMoreThanOneAncestors), ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances2MeanValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
\t\t\twaveFrontDistances2MedianValue = matrix(nrow=1 ,ncol=(wavefrontDistanceSlices+1))
\t\t\tfor (i in 0:wavefrontDistanceSlices)
\t\t\t\t{
\t\t\t\t\ttime = minStartYear+(i*wavefrontDistanceTimeInterval)
\t\t\t\t\tslicedTimes[1,i+1] = time; n = 0
\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tn = n+1
\t\t\t\t\t\t\t\t\twaveFrontDistances1Values[n,i+1] = waveFrontDistances1List[[t]][i+1,2]
\t\t\t\t\t\t\t\t\twaveFrontDistances2Values[n,i+1] = waveFrontDistances2List[[t]][i+1,2]
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tquantiles = quantile(waveFrontDistances1Values[,i], probs=c(0.025,0.975))
\t\t\t\t\tlower_l_1[1,i+1] = as.numeric(quantiles[1])
\t\t\t\t\tupper_l_1[1,i+1] = as.numeric(quantiles[2])
\t\t\t\t\tquantiles = quantile(waveFrontDistances2Values[,i], probs=c(0.025,0.975))
\t\t\t\t\tlower_l_2[1,i+1] = as.numeric(quantiles[1])
\t\t\t\t\tupper_l_2[1,i+1] = as.numeric(quantiles[2])
\t\t\t\t\twaveFrontDistances1MeanValue[1,i+1] = mean(waveFrontDistances1Values[,i+1])
\t\t\t\t\twaveFrontDistances1MedianValue[1,i+1] = median(waveFrontDistances1Values[,i+1])
\t\t\t\t\twaveFrontDistances2MeanValue[1,i+1] = mean(waveFrontDistances2Values[,i+1])
\t\t\t\t\twaveFrontDistances2MedianValue[1,i+1] = median(waveFrontDistances2Values[,i+1])
\t\t\t\t}
\t\t\tyLim1 = c(min(waveFrontDistances1List[[1]][,2]),maxDistance1)
\t\t\tyLim2 = c(min(waveFrontDistances1List[[1]][,2]),maxDistance2)
\t\t\ttreeIDs = paste("distance_tree", treeIDs, sep="")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=2)
\t\t\ttab[,1] = slicedTimes; tab[,2] = waveFrontDistances1MeanValue; colnames(tab) = c("time","distance")\t
\t\t\twrite.table(tab, file=paste(outputName,"_mean_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = waveFrontDistances1MedianValue; colnames(tab) = c("time","distance")\t
\t\t\twrite.table(tab, file=paste(outputName,"_median_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = waveFrontDistances2MeanValue; colnames(tab) = c("time","distance")\t
\t\t\twrite.table(tab, file=paste(outputName,"_mean_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = waveFrontDistances2MedianValue; colnames(tab) = c("time","distance")\t
\t\t\twrite.table(tab, file=paste(outputName,"_median_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = lower_l_1; tab[,3] = upper_l_1; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")\t
\t\t\twrite.table(tab, file=paste(outputName,"_95%HPD_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = lower_l_2; tab[,3] = upper_l_2; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")
\t\t\twrite.table(tab, file=paste(outputName,"_95%HPD_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=(1+dim(waveFrontDistances1Values)[1]))
\t\t\tselectedTreeIDs = treeIDs[which(!seq(1,nberOfExtractionFiles,1)%in%extractionsWithMoreThanOneAncestors)]
\t\t\ttab[,1] = slicedTimes; tab[,2:(1+dim(waveFrontDistances1Values)[1])] = t(waveFrontDistances1Values); colnames(tab) = c("time",selectedTreeIDs)\t
\t\t\twrite.table(tab, file=paste(outputName,"_spatial_wavefront_distances.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2:(1+dim(waveFrontDistances2Values)[1])] = t(waveFrontDistances2Values); colnames(tab) = c("time",selectedTreeIDs)\t
\t\t\twrite.table(tab, file=paste(outputName,"_patristic_wavefront_distances.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\t
\t\t\txLab="time"; yLab = "wavefront distance (km)"
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_spatial_wavefront_distance_1.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Furthest extent of epidemic wavefront (spatial distance from epidemic origin)"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(waveFrontDistances1List[[1]][,1], waveFrontDistances1List[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim1, xlim=xLim)
\t\t\tif (nberOfExtractionFiles > 1)
\t\t\t\t{
\t\t\t\t\tfor (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors) lines(waveFrontDistances1List[[t]][,1], waveFrontDistances1List[[t]][,2], lwd=0.05)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_spatial_wavefront_distance_1.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5,height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_patristic_wavefront_distance_1.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Furthest extent of epidemic wavefront (patristic distance from epidemic origin)"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(waveFrontDistances2List[[1]][,1], waveFrontDistances2List[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim2, xlim=xLim)
\t\t\tif (nberOfExtractionFiles > 1)
\t\t\t\t{
\t\t\t\t\tfor (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (!t%in%extractionsWithMoreThanOneAncestors) lines(waveFrontDistances2List[[t]][,1], waveFrontDistances2List[[t]][,2], lwd=0.05)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_patristic_wavefront_distance_1.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_spatial_wavefront_distance_2.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Furthest extent of epidemic wavefront (spatial distance from epidemic origin)"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(slicedTimes, waveFrontDistances1MedianValue, type="l", axes=F, ann=F, ylim=yLim1, xlim=xLim)
\t\t\txx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
\t\t\tgetOption("scipen"); opt = options("scipen"=20)
\t\t\tpolygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
\t\t\tlines(slicedTimes, waveFrontDistances1MedianValue, lwd=1)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_spatial_wavefront_distance_2.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_patristic_wavefront_distance_2.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of epidemic wavefront (patristic distance from epidemic origin)"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(slicedTimes, waveFrontDistances2MedianValue, type="l", axes=F, ann=F, ylim=yLim2, xlim=xLim)
\t\t\txx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_2,rev(upper_l_2))
\t\t\tgetOption("scipen"); opt = options("scipen"=20)
\t\t\tpolygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
\t\t\tlines(slicedTimes, waveFrontDistances2MedianValue, lwd=1)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_patristic_wavefront_distance_2.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t}
\tif ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&(dispersalVelocityGraph == TRUE))
\t\t{
\t\t\tcat("Building branch dispersal velocity evolution graphs", "\\n", sep="")
\t\t\tslicedTimes = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tlower_l = matrix(nrow=1, ncol=dispersalVelocitySlices); upper_l = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tnumberOfBranchesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
\t\t\tnumberOfBranchesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tnumberOfBranchesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tmeanDispersalVelocitiesMedianValue = matrix(nrow=1 ,ncol=dispersalVelocitySlices)
\t\t\tfor (i in 1:dispersalVelocitySlices)
\t\t\t\t{
\t\t\t\t\tslicedTimes[1,i] = meanDispersalVelocityList[[1]][i,1]
\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tnumberOfBranchesValues[t,i] = numberOfBranchesList[[t]][i,2]
\t\t\t\t\t\t\tmeanDispersalVelocitiesValues[t,i] = meanDispersalVelocityList[[t]][i,2]
\t\t\t\t\t\t}
\t\t\t\t\tnumberOfBranchesMeanValue[1,i] = mean(numberOfBranchesValues[,i], na.rm=T)
\t\t\t\t\tnumberOfBranchesMedianValue[1,i] = median(numberOfBranchesValues[,i], na.rm=T)
\t\t\t\t\tquantiles = quantile(meanDispersalVelocitiesValues[,i], probs=c(0.025,0.975), na.rm=T)
\t\t\t\t\tlower_l[1,i] = as.numeric(quantiles[1]); upper_l[1,i] = as.numeric(quantiles[2])
\t\t\t\t\tmeanDispersalVelocitiesMeanValue[1,i] = mean(meanDispersalVelocitiesValues[,i], na.rm=T)
\t\t\t\t\tmeanDispersalVelocitiesMedianValue[1,i] = median(meanDispersalVelocitiesValues[,i], na.rm=T)
\t\t\t\t}
\t\t\tyLim = c(0, max(upper_l, na.rm=T))
\t\t\ttreeIDs = paste("distance_tree", treeIDs, sep="")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = meanDispersalVelocitiesMeanValue
\t\t\ttab[,3] = numberOfBranchesMeanValue; colnames(tab) = c("time","velocity","number_of_branches")\t
\t\t\twrite.table(tab, file=paste(outputName,"_mean_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = meanDispersalVelocitiesMedianValue
\t\t\ttab[,3] = numberOfBranchesMedianValue; colnames(tab) = c("time","velocity","number_of_branches")\t
\t\t\twrite.table(tab, file=paste(outputName,"_median_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")\t
\t\t\twrite.table(tab, file=paste(outputName,"_95%HPD_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")

\t\t\txLab = "time"; yLab = "weighted branch dispersal velocity"
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_mean_branch_dispersal_velocity_1.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of mean branch dispersal velocity"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(meanDispersalVelocityList[[1]][,1], meanDispersalVelocityList[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\tif (nberOfExtractionFiles > 1)
\t\t\t\t{
\t\t\t\t\tfor (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tlines(meanDispersalVelocityList[[t]][,1],meanDispersalVelocityList[[t]][,2],lwd=LWD)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_1.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t\t\t\t
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_mean_branch_dispersal_velocity_2.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of mean branch dispersal velocity"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(slicedTimes, meanDispersalVelocitiesMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\txx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
\t\t\tgetOption("scipen"); opt = options("scipen"=20)
\t\t\tpolygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
\t\t\tlines(slicedTimes, meanDispersalVelocitiesMedianValue, lwd=1)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_2.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tslicedTimes = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tlower_l = matrix(nrow=1, ncol=dispersalVelocitySlices); upper_l = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
\t\t\tweightedDispersalVelocitiesMedianValue = matrix(nrow=1 ,ncol=dispersalVelocitySlices)
\t\t\tfor (i in 1:dispersalVelocitySlices)
\t\t\t\t{
\t\t\t\t\tslicedTimes[1,i] = weightedDispersalVelocityList[[1]][i,1]
\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tweightedDispersalVelocitiesValues[t,i] = weightedDispersalVelocityList[[t]][i,2]
\t\t\t\t\t\t}
\t\t\t\t\tquantiles = quantile(weightedDispersalVelocitiesValues[,i], probs=c(0.025,0.975), na.rm=T)
\t\t\t\t\tlower_l[1,i] = as.numeric(quantiles[1]); upper_l[1,i] = as.numeric(quantiles[2])
\t\t\t\t\tweightedDispersalVelocitiesMeanValue[1,i] = mean(weightedDispersalVelocitiesValues[,i], na.rm=T)
\t\t\t\t\tweightedDispersalVelocitiesMedianValue[1,i] = median(weightedDispersalVelocitiesValues[,i], na.rm=T)
\t\t\t\t}
\t\t\tyLim = c(0, max(upper_l, na.rm=T))
\t\t\ttreeIDs = paste("distance_tree", treeIDs, sep="")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = weightedDispersalVelocitiesMeanValue
\t\t\ttab[,3] = numberOfBranchesMeanValue; colnames(tab) = c("time","velocity","number_of_branches")\t
\t\t\twrite.table(tab, file=paste(outputName,"_mean_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab[,1] = slicedTimes; tab[,2] = weightedDispersalVelocitiesMedianValue
\t\t\ttab[,3] = numberOfBranchesMedianValue; colnames(tab) = c("time","velocity","number_of_branches")\t
\t\t\twrite.table(tab, file=paste(outputName,"_median_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")\t
\t\t\twrite.table(tab, file=paste(outputName,"_95%HPD_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\\t")

\t\t\txLab = "time"; yLab = "weighted branch dispersal velocity"
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_1.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of weighted branch dispersal velocity"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(weightedDispersalVelocityList[[1]][,1], weightedDispersalVelocityList[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\tif (nberOfExtractionFiles > 1)
\t\t\t\t{
\t\t\t\t\tfor (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tlines(weightedDispersalVelocityList[[t]][,1],weightedDispersalVelocityList[[t]][,2],lwd=LWD)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_1.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\t
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_2.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of weighted branch dispersal velocity"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(slicedTimes, weightedDispersalVelocitiesMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\txx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
\t\t\tgetOption("scipen"); opt = options("scipen"=20)
\t\t\tpolygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
\t\t\tlines(slicedTimes, weightedDispersalVelocitiesMedianValue, lwd=1)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_2.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t}
}
