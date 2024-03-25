spreadStatistics = function(localTreesDirectory="", nberOfExtractionFiles=1, timeSlices=200, onlyTipBranches=F, showingPlots=TRUE, outputName=gsub(" ","_",date()), nberOfCores=1, slidingWindow=NA, simulations=FALSE, discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot=FALSE) {

	nberOfStatistics = 9; treeIDs = c()
	registerDoMC(cores=nberOfCores)
  		# (1) mean branch dispersal velocity
  		# (2) weighted branch dispersal velocity
  		# (3) branch dispersal velocity variation among branches (CV)
  		# (4) original diffuson coefficient (Pybus et al. 2012)
  		# (5) weighted diffuson coefficient (Trovao et al. 2015)
  		# (6) diffusion coefficient variation among branches (CV)
  		# (7) isolation-by-distance (IBD) signal estimated by r_S
  		#	  (Spearman correlation between patristic and geographic distances)
  		# (8) isolation-by-distance (IBD) signal estimated by r_P
  		#	  (Pearson correlation between patristic and geographic distances)
  		# (9) isolation-by-distance (IBD) signal estimated by r_P
  		#	  (Pearson correlation between patristic and log-transformed geographic distances)
	meanStatistics = matrix(nrow=(nberOfExtractionFiles), ncol=nberOfStatistics)
	branchVelocities = c() # not used, just to obtain an overall distribution of velocities
	sd_var_velocity = matrix(nrow=(nberOfExtractionFiles), ncol=2) # not used either
	medianMeanStatistics = matrix(nrow=1, ncol=nberOfStatistics)
	ciMeanStatistics = matrix(nrow=2, ncol=nberOfStatistics)
	waveFrontDistances1List = list()
	waveFrontDistances2List = list()
	meanBranchDispersalVelocityList = list()
	weightedBranchDispersalVelocityList = list()
	numberOfBranchesList = list()
	dispersalOrientationList = list()
	cat("Estimation of summary statistics", "\n", sep="")
	onlyOneAncestor = TRUE; extractionsWithMoreThanOneAncestors = c()
	if ((timeSlices==0)|is.na(timeSlices)|is.null(timeSlices)) onlyOneAncestor = FALSE
	dispersalVelocityGraph = FALSE
	if (!is.na(slidingWindow)) dispersalVelocityGraph = TRUE
	for (t in 1:nberOfExtractionFiles)
		{
			if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
			if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear,startYear)),]
			if (sum(!data[,"node1"]%in%data[,"node2"]) > 2)
				{
					onlyOneAncestor = FALSE; extractionsWithMoreThanOneAncestors = c(extractionsWithMoreThanOneAncestors, t)
				}
			treeIDs = c(treeIDs, data[1,"treeID"])
			if (onlyTipBranches == TRUE)
				{
					indices = c()
					for (i in 1:dim(data)[1])
						{
							n = data[i,"node2"]
							if (length(data[data[,"node1"]==n,"node1"]) == 0)
								{
									indices = c(indices, i)
								}
						}
					data = data[indices,]	
				}
			nberOfConnections = dim(data)[1]
			if (t == 1)
				{
					minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					minStartYear = min(data[,"startYear"])
		    		minEndYear = min(data[,"startYear"])
		    		maxEndYear = max(data[,"endYear"])
				}	else	{
					if (minLon > min(min(data[,"endLon"]),min(data[,"startLon"]))) minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					if (maxLon < max(max(data[,"endLon"]),max(data[,"startLon"]))) maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					if (minLat > min(min(data[,"endLat"]),min(data[,"startLat"]))) minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					if (maxLat < max(max(data[,"endLat"]),max(data[,"startLat"]))) maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					if (minStartYear > min(data[,"startYear"])) minStartYear = min(data[,"startYear"])
					if (minEndYear > min(data[,"endYear"])) minEndYear = min(data[,"endYear"])
					if (maxEndYear < max(data[,"endYear"])) maxEndYear = max(data[,"endYear"])
				}
		}
	if ((onlyOneAncestor == FALSE)&(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE))
		{
			cat("Discarding",length(extractionsWithMoreThanOneAncestors),"extraction tables with more without a single","\n")
			cat("\t","common ancestor for generating wavefront plots","\n")
		}
	wavefrontDistanceSlices = timeSlices
	wavefrontDistanceTimeInterval = (maxEndYear-minStartYear)/wavefrontDistanceSlices
	dispersalVelocitySlices = timeSlices
	dispersalVelocityTimeInterval = (maxEndYear-minStartYear)/dispersalVelocitySlices
	xLim = c(minStartYear, maxEndYear)
	for (t in 1:nberOfExtractionFiles)
		{
			if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
			if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear,startYear)),]
			nberOfConnections = dim(data)[1]
			distances = matrix(0, nrow=dim(data)[1], ncol=1)
	 		colnames(distances) = c("greatCircleDist_km")
			for (i in 1:nberOfConnections)
				{
					x1 = cbind(data[i,"startLon"],data[i,"startLat"])
	 				x2 = cbind(data[i,"endLon"],data[i,"endLat"])
	 				distances[i] = rdist.earth(x1, x2, miles=F, R=NULL)
				}
			if (length(grep("greatCircleDist_km",colnames(data))) > 0)
				{
					data[,"greatCircleDist_km"] = distances
				}	else	{
					data = cbind(data,distances)
				}
			branchMeasures = matrix(nrow=nberOfConnections, ncol=2)
			weightedBranchDispersalVelocity_numerator = 0; weightedBranchDispersalVelocity_denominator = 0
			weightedDiffusionCoefficient_numerator = 0; weightedDiffusionCoefficient_denominator = 0
			for (i in 1:nberOfConnections)
				{
					dispersalTime = data[i,"endYear"]-data[i,"startYear"]
		    		BranchDispersalVelocity = data[i,"greatCircleDist_km"]/(dispersalTime)
		    		branchOriginalDiffusionCoefficient = (data[i,"greatCircleDist_km"]^2)/(4*dispersalTime)
		    		branchMeasures[i,1] = BranchDispersalVelocity
		    		branchMeasures[i,2] = branchOriginalDiffusionCoefficient
		    		weightedBranchDispersalVelocity_numerator = weightedBranchDispersalVelocity_numerator + data[i,"greatCircleDist_km"]
		    		weightedBranchDispersalVelocity_denominator = weightedBranchDispersalVelocity_denominator + dispersalTime
		    		weightedDiffusionCoefficient_numerator = weightedDiffusionCoefficient_numerator + (data[i,"greatCircleDist_km"]^2)
		    		weightedDiffusionCoefficient_denominator = weightedDiffusionCoefficient_denominator + (4*dispersalTime)
		  		}
		  	branchVelocities = c(branchVelocities, branchMeasures[,1])
			sd_var_velocity[t,] = cbind(sd(branchMeasures[,1]), var(branchMeasures[,1]))
			meanStatistics[t,1] = mean(branchMeasures[,1])
			meanStatistics[t,2] = weightedBranchDispersalVelocity_numerator/weightedBranchDispersalVelocity_denominator
			meanStatistics[t,3] = sd(branchMeasures[,1])/mean(branchMeasures[,1])
			meanStatistics[t,4] = mean(branchMeasures[,2])
		    meanStatistics[t,5] = weightedDiffusionCoefficient_numerator/weightedDiffusionCoefficient_denominator
			meanStatistics[t,6] = sd(branchMeasures[,2])/mean(branchMeasures[,2])
			tipNodeIndices = which(!data[,"node2"]%in%data[,"node1"])
			distTree = matrix(nrow=length(tipNodeIndices), ncol=length(tipNodeIndices))
			for (k in 2:dim(distTree)[1])
				{
					for (l in 1:(k-1))
						{
							index1 = tipNodeIndices[k]
							index2 = tipNodeIndices[l]
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
							distTree[k,l] = patristic_dis; distTree[l,k] = patristic_dis
						}
				}
			distsGeo = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
			for (k in 2:dim(distsGeo)[1])
				{
					for (l in 1:(k-1))
						{
							index1 = tipNodeIndices[k]
							index2 = tipNodeIndices[l]
							x1 = cbind(data[index1,"endLon"], data[index1,"endLat"])
							x2 = cbind(data[index2,"endLon"], data[index2,"endLat"])
							distsGeo[k,l] = rdist.earth(x1, x2, miles=F, R=NULL)
							distsGeo[l,k] = distsGeo[k,l]
						}
				}			
			meanStatistics[t,7] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="spearman") # r_S
			meanStatistics[t,8] = cor(distTree[lower.tri(distTree)],distsGeo[lower.tri(distsGeo)], method="pearson") # r_P #1
			meanStatistics[t,9] = cor(distTree[lower.tri(distTree)],log(distsGeo[lower.tri(distsGeo)]+1), method="pearson") # r_P #2
		}
	if ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&((onlyOneAncestor == TRUE)|(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE)))
		{
			buffer = list()
			buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {				
			# for (t in 1:nberOfExtractionFiles) {
					if (!t%in%extractionsWithMoreThanOneAncestors)
						{
							if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
							if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
							data = read.csv(fileName, h=T)
							data = data[with(data, order(endYear, startYear)),]
							nberOfConnections = dim(data)[1]
							waveFrontDistances1 = matrix(nrow=wavefrontDistanceSlices+1, ncol=2) # spatial wavefront distance
							waveFrontDistances2 = matrix(nrow=wavefrontDistanceSlices+1, ncol=2) # patristic wavefront distance
							waveFrontDistances1[,1] = seq(minStartYear,maxEndYear,wavefrontDistanceTimeInterval)
							waveFrontDistances2[,1] = seq(minStartYear,maxEndYear,wavefrontDistanceTimeInterval)
							data = data[order(data[,"endYear"]),]
							ancestralIndex = which(data[,"startYear"]==min(data[,"startYear"]))[1]
							ancestralNode = data[ancestralIndex,"node1"]
							startingYear = data[ancestralIndex,"startYear"]
							originLocationX = data[ancestralIndex,"startLon"]
							originLocationY = data[ancestralIndex,"startLat"]
							originLocation = cbind(originLocationX, originLocationY)
							maxDistance1 = 0; maxDistance2 = 0
							for (i in 0:wavefrontDistanceSlices)
								{
									time = minStartYear+(i*wavefrontDistanceTimeInterval)
									for (j in 1:nberOfConnections)
										{
											if ((data[j,"startYear"]<time)&(time<data[j,"endYear"]))
												{
													timeProportion = (time-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
													pointLocationX = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion)
													pointLocationY = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion)
													pointLocation = cbind(pointLocationX, pointLocationY)
													distance1 = rdist.earth(originLocation, pointLocation, miles=F, R=NULL)				    		 		
						    							coordinatesOfStartNode = cbind(data[j,"startLon"], data[j,"startLat"])
						    							distance2 = rdist.earth(coordinatesOfStartNode, pointLocation, miles=F, R=NULL)
						    							originNode = data[j,"node1"]
						    							while (originNode != ancestralNode)
						    								{
						    									index = which(data[,"node2"]==originNode)[1]
						    									coordinatesOfEndNode = cbind(data[index,"endLon"], data[index,"endLat"])
						    									coordinatesOfStartNode = cbind(data[index,"startLon"], data[index,"startLat"])
						    									distance2 = distance2 + rdist.earth(coordinatesOfStartNode, coordinatesOfEndNode, miles=F, R=NULL)
						    									originNode = data[index,"node1"]
						    								}
													if (distance1 > maxDistance1) maxDistance1 = distance1
						    							if (distance2 > maxDistance2) maxDistance2 = distance2
												}
											
										}
									waveFrontDistances1[i+1,2] = maxDistance1
						    		waveFrontDistances2[i+1,2] = maxDistance2
								}
							colnames(waveFrontDistances1) = c("endYear", "waveFrontDistances1")
							colnames(waveFrontDistances2) = c("endYear", "waveFrontDistances2")
							waveFrontDistances = list()
							waveFrontDistances[[1]] = waveFrontDistances1
							waveFrontDistances[[2]] = waveFrontDistances2
							# buffer[[t]] = waveFrontDistances
							waveFrontDistances
						}
				}
			for (t in 1:length(buffer))
				{
					if (!t%in%extractionsWithMoreThanOneAncestors)
						{
							waveFrontDistances1List[[t]] = buffer[[t]][[1]]
							waveFrontDistances2List[[t]] = buffer[[t]][[2]]
						}
				}
			maxDistance1 = 0; maxDistance2 = 0
			for (t in 1:nberOfExtractionFiles)
				{
					if (!t%in%extractionsWithMoreThanOneAncestors)
						{
							if (max(waveFrontDistances1List[[t]][,2]) > maxDistance1) maxDistance1 = max(waveFrontDistances1List[[t]][,2])
							if (max(waveFrontDistances2List[[t]][,2]) > maxDistance2) maxDistance2 = max(waveFrontDistances2List[[t]][,2])
						}
				}
		}
	if ((onlyTipBranches == FALSE)&(dispersalVelocityGraph == TRUE)&(nberOfExtractionFiles > 1))
		{
			startEndTimes = matrix(nrow=dispersalVelocitySlices, ncol=3)
			for (i in 1:dispersalVelocitySlices)
				{
					time = minStartYear+((i-1)*dispersalVelocityTimeInterval)+(dispersalVelocityTimeInterval/2)
					startTime = time - (slidingWindow/2)
					endTime = time + (slidingWindow/2)
					startEndTimes[i,1:3] = cbind(time, startTime, endTime)
				}
			buffer = list()
			buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {		
			# for (t in 1:nberOfExtractionFiles) {
					if (simulations == FALSE) fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
					if (simulations == TRUE) fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
					data = read.csv(fileName, h=T)
					data = data[with(data, order(endYear, startYear)),]
					nberOfConnections = dim(data)[1]
					meanBranchDispersalVelocities = matrix(nrow=dispersalVelocitySlices, ncol=2)
					weightedBranchDispersalVelocities = matrix(nrow=dispersalVelocitySlices, ncol=2)
					numberOfBranches = matrix(nrow=dispersalVelocitySlices, ncol=2)
					data = data[order(data[,"endYear"]),]
					for (i in 1:dispersalVelocitySlices)
						{
							n = 0; vS = 0; dS = 0; tS = 0 
							time = startEndTimes[i,1]
							startTime = startEndTimes[i,2]
							endTime = startEndTimes[i,3]
							for (j in 1:nberOfConnections)
								{
									branchInwavefrontDistanceTimeInterval = FALSE
									if ((data[j,"startYear"]<startTime)&(data[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
									if ((data[j,"startYear"]>startTime)&(data[j,"startYear"]<endTime)) branchInwavefrontDistanceTimeInterval = TRUE
									if ((data[j,"endYear"]>startTime)&(data[j,"endYear"]<endTime)) branchInwavefrontDistanceTimeInterval = TRUE
									if (branchInwavefrontDistanceTimeInterval == TRUE)
										{
											if (data[j,"startYear"] > startTime) time1 = data[j,"startYear"]
											if (data[j,"startYear"] <= startTime) time1 = startTime
											if (data[j,"endYear"] < endTime) time2 = data[j,"endYear"]
											if (data[j,"endYear"] >= endTime) time2 = endTime
											branchTimeInInterval = time2 - time1
											timeProportion1 = (time1-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
											pointLocation1X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion1)
											pointLocation1Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion1)
											pointLocation1 = cbind(pointLocation1X, pointLocation1Y)
											timeProportion2 = (time2-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
											pointLocation2X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion2)
											pointLocation2Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion2)
											pointLocation2 = cbind(pointLocation2X, pointLocation2Y)
											branchDistInInterval = rdist.earth(pointLocation1, pointLocation2, miles=F, R=NULL)				    		 		
				    							n = n+1; vS = vS + (branchDistInInterval/branchTimeInInterval)
				    							dS = dS + branchDistInInterval; tS = tS + branchTimeInInterval
										}
								}
							meanBranchDispersalVelocities[i,1] = time
							meanBranchDispersalVelocities[i,2] = vS/n
							weightedBranchDispersalVelocities[i,1] = time
							weightedBranchDispersalVelocities[i,2] = dS/tS
							numberOfBranches[i,1] = time
							numberOfBranches[i,2] = n
						}
					colnames(meanBranchDispersalVelocities) = c("year","meanBranchDispersalVelocity")
					colnames(weightedBranchDispersalVelocities) = c("year","weightedBranchDispersalVelocity")
					dispersalVelocities = list()
					dispersalVelocities[[1]] = meanBranchDispersalVelocities
					dispersalVelocities[[2]] = weightedBranchDispersalVelocities
					dispersalVelocities[[3]] = numberOfBranches
					# buffer[[t]] = dispersalVelocities
					dispersalVelocities
				}
			for (t in 1:length(buffer))
				{
					meanBranchDispersalVelocityList[[t]] = buffer[[t]][[1]]
					weightedBranchDispersalVelocityList[[t]] = buffer[[t]][[2]]
					numberOfBranchesList[[t]] = buffer[[t]][[3]]
				}
		}
	for (i in 1:nberOfStatistics)
		{
			medianMeanStatistics[1,i] = median(meanStatistics[,i], na.rm=T)
			quantiles = quantile(meanStatistics[,i], probs=c(0.025,0.975), na.rm=T)
			ciMeanStatistics[1,i] = as.numeric(quantiles[1]); ciMeanStatistics[2,i] = as.numeric(quantiles[2])
			HPD = HDInterval::hdi(meanStatistics[,i])[1:2]
			ciMeanStatistics[1,i] = as.numeric(HPD[1]); ciMeanStatistics[2,i] = as.numeric(HPD[2])
		}
	cat("Median value of the mean branch dispersal velocity = ",medianMeanStatistics[1,1],"\n	95% HPD = [",ciMeanStatistics[1,1],", ",ciMeanStatistics[2,1],"]","\n",sep="")
	cat("Median value of the weighted branch dispersal velocity = ",medianMeanStatistics[1,2],"\n	95% HPD = [",ciMeanStatistics[1,2],", ",ciMeanStatistics[2,2],"]","\n",sep="")	
	cat("Median value of the original diffusion coefficient = ",medianMeanStatistics[1,4],"\n	95% HPD = [",ciMeanStatistics[1,4],", ",ciMeanStatistics[2,4],"]","\n",sep="")	
	cat("Median value of the weighted diffusion coefficient = ",medianMeanStatistics[1,5],"\n	95% HPD = [",ciMeanStatistics[1,5],", ",ciMeanStatistics[2,5],"]","\n",sep="")	
	cat("Median value of the isolation-by-distance (IBD) signal (rS) = ",medianMeanStatistics[1,7],"\n	95% HPD = [",ciMeanStatistics[1,7],", ",ciMeanStatistics[2,7],"]","\n",sep="")	
	cat("Median value of the isolation-by-distance (IBD) signal (rP #1) = ",medianMeanStatistics[1,8],"\n	95% HPD = [",ciMeanStatistics[1,8],", ",ciMeanStatistics[2,8],"]","\n",sep="")	
	cat("Median value of the isolation-by-distance (IBD) signal (rP #2) = ",medianMeanStatistics[1,9],"\n	95% HPD = [",ciMeanStatistics[1,9],", ",ciMeanStatistics[2,9],"]","\n",sep="")	
	colnames(meanStatistics) = c("mean_branch_dispersal_velocity", "weighted_branch_dispersal_velocity", "branch_dispersal_velocity_variation_among_branches",
								 "original_diffusion_coefficient", "weighted_diffusion_coefficient", "diffusion_coefficient_variation_among_branches",
								 "isolation_by_distance_signal_rS", "isolation_by_distance_signal_rP1", "isolation_by_distance_signal_rP2")
	write.table(meanStatistics, file=paste(outputName,"_estimated_dispersal_statistics.txt",sep=""), quote=F, row.names=F, sep="\t")

	LWD = 0.2
	if (nberOfExtractionFiles > 1)
		{
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName, "_mean_branch_dispersal_velocity_variation.pdf",sep=""), width=5, height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,1]), meanStatistics[,3]))
			# kde = kde(cbind(log10(meanStatistics[,1]), meanStatistics[,3]), H=H)
			H = Hpi(cbind(meanStatistics[,1], meanStatistics[,3]))
			kde = kde(cbind(meanStatistics[,1], meanStatistics[,3]), H=H)
			text = "Kernel density estimates of mean branch dispersal velocity parameters"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean branch velocity"; yLab = "mean branch velocity variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()

			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""),width=5,height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,2]), meanStatistics[,3]))
			# kde = kde(cbind(log10(meanStatistics[,2]), meanStatistics[,3]), H=H)
			H = Hpi(cbind(meanStatistics[,2], meanStatistics[,3]))
			kde = kde(cbind(meanStatistics[,2], meanStatistics[,3]), H=H)
			text = "Kernel density estimates of weighted branch dispersal velocity parameters"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "weighted dispersal velocity"; yLab="weighted dispersal velocity variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()

			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""),width=5,height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,4]), meanStatistics[,6]))
			# kde = kde(cbind(log10(meanStatistics[,4]), meanStatistics[,6]), H=H)
			H = Hpi(cbind(meanStatistics[,4], meanStatistics[,6]))
			kde = kde(cbind(meanStatistics[,4], meanStatistics[,6]), H=H)
			text1 = "Kernel density estimates of original diffusion coefficient parameters"; text2 = "(Pybus et al. 2012)"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean original diffusion coefficient"; yLab="diffusion coefficient variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text1, cex.main=0.6, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_original_diffusion_coefficient_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""), width=5, height=5)
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0), lwd=LWD, col="gray30")
			# H = Hpi(cbind(log10(meanStatistics[,5]), meanStatistics[,6]))
			# kde = kde(cbind(log10(meanStatistics[,5]), meanStatistics[,6]), H=H)
			H = Hpi(cbind(meanStatistics[,5],meanStatistics[,6]))
			kde = kde(cbind(meanStatistics[,5],meanStatistics[,6]),H=H)
			text1 = "Kernel density estimates of weighted diffusion coefficient parameters"; text2 = "(Trovao et al. 2015)"
			colours = c("#FFFFFF","#D2D3D3","#9D9FA3","#6A6A6D")
			xLab = "mean weighted diffusion coefficient"; yLab="diffusion coefficient variation among branches"
			plot(kde, display="filled.contour2", cont=c(50,75,95), col=colours, axes=F, ann=F)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text1, cex.main=0.6, col.main="gray30")
			legend("bottomright", c("","95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=15, pt.cex=1.5, col=colours, box.lty=0, cex=0.6, y.intersp=1.2)
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_diffusion_coefficient_variation.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
		}

	if ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&((onlyOneAncestor == TRUE)|(discardExtractionTablesWithMoreThanOneAncestorForWavefrontPlot == TRUE)))
		{
			cat("Building wavefront distance evolution graphs", "\n", sep="")
			wavefrontDistanceTimeInterval = (maxEndYear-minStartYear)/wavefrontDistanceSlices
			slicedTimes = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			lower_l_1 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			upper_l_1 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			lower_l_2 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			upper_l_2 = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances1Values = matrix(nrow=nberOfExtractionFiles-length(extractionsWithMoreThanOneAncestors), ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances1MeanValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances1MedianValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances2Values = matrix(nrow=nberOfExtractionFiles-length(extractionsWithMoreThanOneAncestors), ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances2MeanValue = matrix(nrow=1, ncol=(wavefrontDistanceSlices+1))
			waveFrontDistances2MedianValue = matrix(nrow=1 ,ncol=(wavefrontDistanceSlices+1))
			for (i in 0:wavefrontDistanceSlices)
				{
					time = minStartYear+(i*wavefrontDistanceTimeInterval)
					slicedTimes[1,i+1] = time; n = 0
					for (t in 1:nberOfExtractionFiles)
						{
							if (!t%in%extractionsWithMoreThanOneAncestors)
								{
									n = n+1
									waveFrontDistances1Values[n,i+1] = waveFrontDistances1List[[t]][i+1,2]
									waveFrontDistances2Values[n,i+1] = waveFrontDistances2List[[t]][i+1,2]
								}
						}
					quantiles = quantile(waveFrontDistances1Values[,i], probs=c(0.025,0.975))
					lower_l_1[1,i+1] = as.numeric(quantiles[1]); upper_l_1[1,i+1] = as.numeric(quantiles[2])
					HPD = HDInterval::hdi(waveFrontDistances1Values[,i])[1:2]
					lower_l_1[1,i+1] = as.numeric(HPD[1]); upper_l_1[1,i+1] = as.numeric(HPD[2])
					quantiles = quantile(waveFrontDistances2Values[,i], probs=c(0.025,0.975))
					lower_l_2[1,i+1] = as.numeric(quantiles[1]); upper_l_2[1,i+1] = as.numeric(quantiles[2])
					HPD = HDInterval::hdi(waveFrontDistances2Values[,i])[1:2]
					lower_l_2[1,i+1] = as.numeric(HPD[1]); upper_l_2[1,i+1] = as.numeric(HPD[2])
					waveFrontDistances1MeanValue[1,i+1] = mean(waveFrontDistances1Values[,i+1])
					waveFrontDistances1MedianValue[1,i+1] = median(waveFrontDistances1Values[,i+1])
					waveFrontDistances2MeanValue[1,i+1] = mean(waveFrontDistances2Values[,i+1])
					waveFrontDistances2MedianValue[1,i+1] = median(waveFrontDistances2Values[,i+1])
				}
			yLim1 = c(min(waveFrontDistances1List[[1]][,2]),maxDistance1)
			yLim2 = c(min(waveFrontDistances1List[[1]][,2]),maxDistance2)
			treeIDs = paste("distance_tree", treeIDs, sep="")
			tab = matrix(nrow=length(slicedTimes), ncol=2)
			tab[,1] = slicedTimes; tab[,2] = waveFrontDistances1MeanValue; colnames(tab) = c("time","distance")	
			write.table(tab, file=paste(outputName,"_mean_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = waveFrontDistances1MedianValue; colnames(tab) = c("time","distance")	
			write.table(tab, file=paste(outputName,"_median_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = waveFrontDistances2MeanValue; colnames(tab) = c("time","distance")	
			write.table(tab, file=paste(outputName,"_mean_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = waveFrontDistances2MedianValue; colnames(tab) = c("time","distance")	
			write.table(tab, file=paste(outputName,"_median_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = lower_l_1; tab[,3] = upper_l_1; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")	
			write.table(tab, file=paste(outputName,"_95%HPD_spatial_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = lower_l_2; tab[,3] = upper_l_2; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")
			write.table(tab, file=paste(outputName,"_95%HPD_patristic_wavefront_distance.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab = matrix(nrow=length(slicedTimes), ncol=(1+dim(waveFrontDistances1Values)[1]))
			selectedTreeIDs = treeIDs[which(!seq(1,nberOfExtractionFiles,1)%in%extractionsWithMoreThanOneAncestors)]
			tab[,1] = slicedTimes; tab[,2:(1+dim(waveFrontDistances1Values)[1])] = t(waveFrontDistances1Values); colnames(tab) = c("time",selectedTreeIDs)	
			write.table(tab, file=paste(outputName,"_spatial_wavefront_distances.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2:(1+dim(waveFrontDistances2Values)[1])] = t(waveFrontDistances2Values); colnames(tab) = c("time",selectedTreeIDs)	
			write.table(tab, file=paste(outputName,"_patristic_wavefront_distances.txt",sep=""), row.names=F, quote=F, sep="\t")
			
			xLab="time"; yLab = "wavefront distance (km)"
			
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_spatial_wavefront_distance_1.pdf",sep=""), width=5, height=5)
			text = "Furthest extent of epidemic wavefront (spatial distance from epidemic origin)"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(waveFrontDistances1List[[1]][,1], waveFrontDistances1List[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim1, xlim=xLim)
			if (nberOfExtractionFiles > 1)
				{
					for (t in 2:nberOfExtractionFiles)
						{
							if (!t%in%extractionsWithMoreThanOneAncestors) lines(waveFrontDistances1List[[t]][,1], waveFrontDistances1List[[t]][,2], lwd=0.05)
						}
				}
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_spatial_wavefront_distance_1.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5,height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_patristic_wavefront_distance_1.pdf",sep=""), width=5, height=5)
			text = "Furthest extent of epidemic wavefront (patristic distance from epidemic origin)"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(waveFrontDistances2List[[1]][,1], waveFrontDistances2List[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim2, xlim=xLim)
			if (nberOfExtractionFiles > 1)
				{
					for (t in 2:nberOfExtractionFiles)
						{
							if (!t%in%extractionsWithMoreThanOneAncestors) lines(waveFrontDistances2List[[t]][,1], waveFrontDistances2List[[t]][,2], lwd=0.05)
						}
				}
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_patristic_wavefront_distance_1.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_spatial_wavefront_distance_2.pdf",sep=""), width=5, height=5)
			text = "Furthest extent of epidemic wavefront (spatial distance from epidemic origin)"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(slicedTimes, waveFrontDistances1MedianValue, type="l", axes=F, ann=F, ylim=yLim1, xlim=xLim)
			xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
			getOption("scipen"); opt = options("scipen"=20)
			polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
			lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_spatial_wavefront_distance_2.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_patristic_wavefront_distance_2.pdf",sep=""), width=5, height=5)
			text = "Evolution of epidemic wavefront (patristic distance from epidemic origin)"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(slicedTimes, waveFrontDistances2MedianValue, type="l", axes=F, ann=F, ylim=yLim2, xlim=xLim)
			xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_2,rev(upper_l_2))
			getOption("scipen"); opt = options("scipen"=20)
			polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
			lines(slicedTimes, waveFrontDistances2MedianValue, lwd=1)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.55, col.main="gray30"); box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_patristic_wavefront_distance_2.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
		}
	if ((nberOfExtractionFiles > 1)&(onlyTipBranches == FALSE)&(dispersalVelocityGraph == TRUE))
		{
			cat("Building branch dispersal velocity evolution graphs", "\n", sep="")
			slicedTimes = matrix(nrow=1, ncol=dispersalVelocitySlices)
			lower_l = matrix(nrow=1, ncol=dispersalVelocitySlices); upper_l = matrix(nrow=1, ncol=dispersalVelocitySlices)
			numberOfBranchesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
			numberOfBranchesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			numberOfBranchesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			meanBranchDispersalVelocitiesMedianValue = matrix(nrow=1 ,ncol=dispersalVelocitySlices)
			for (i in 1:dispersalVelocitySlices)
				{
					slicedTimes[1,i] = meanBranchDispersalVelocityList[[1]][i,1]
					for (t in 1:nberOfExtractionFiles)
						{
							numberOfBranchesValues[t,i] = numberOfBranchesList[[t]][i,2]
							meanBranchDispersalVelocitiesValues[t,i] = meanBranchDispersalVelocityList[[t]][i,2]
						}
					numberOfBranchesMeanValue[1,i] = mean(numberOfBranchesValues[,i], na.rm=T)
					numberOfBranchesMedianValue[1,i] = median(numberOfBranchesValues[,i], na.rm=T)
					quantiles = quantile(meanBranchDispersalVelocitiesValues[,i], probs=c(0.025,0.975), na.rm=T)
					lower_l[1,i] = as.numeric(quantiles[1]); upper_l[1,i] = as.numeric(quantiles[2])
					HPD = HDInterval::hdi(meanBranchDispersalVelocitiesValues[,i])[1:2]
					lower_l[1,i] = as.numeric(HPD[1]); upper_l[1,i] = as.numeric(HPD[2])
					meanBranchDispersalVelocitiesMeanValue[1,i] = mean(meanBranchDispersalVelocitiesValues[,i], na.rm=T)
					meanBranchDispersalVelocitiesMedianValue[1,i] = median(meanBranchDispersalVelocitiesValues[,i], na.rm=T)
				}
			yLim = c(0, max(upper_l, na.rm=T))
			treeIDs = paste("distance_tree", treeIDs, sep="")
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = meanBranchDispersalVelocitiesMeanValue
			tab[,3] = numberOfBranchesMeanValue; colnames(tab) = c("time","velocity","number_of_branches")	
			write.table(tab, file=paste(outputName,"_mean_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = meanBranchDispersalVelocitiesMedianValue
			tab[,3] = numberOfBranchesMedianValue; colnames(tab) = c("time","velocity","number_of_branches")	
			write.table(tab, file=paste(outputName,"_median_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")	
			write.table(tab, file=paste(outputName,"_95%HPD_mean_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")

			xLab = "time"; yLab = "mean branch dispersal velocity"
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_mean_branch_dispersal_velocity_1.pdf",sep=""), width=5, height=5)
			text = "Evolution of mean branch dispersal velocity"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(meanBranchDispersalVelocityList[[1]][,1], meanBranchDispersalVelocityList[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
			if (nberOfExtractionFiles > 1)
				{
					for (t in 2:nberOfExtractionFiles)
						{
							lines(meanBranchDispersalVelocityList[[t]][,1],meanBranchDispersalVelocityList[[t]][,2],lwd=LWD)
						}
				}
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_1.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
						
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_mean_branch_dispersal_velocity_2.pdf",sep=""), width=5, height=5)
			text = "Evolution of mean branch dispersal velocity"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(slicedTimes, meanBranchDispersalVelocitiesMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
			xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
			getOption("scipen"); opt = options("scipen"=20)
			polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
			lines(slicedTimes, meanBranchDispersalVelocitiesMedianValue, lwd=1)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_mean_branch_dispersal_velocity_2.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			slicedTimes = matrix(nrow=1, ncol=dispersalVelocitySlices)
			lower_l = matrix(nrow=1, ncol=dispersalVelocitySlices); upper_l = matrix(nrow=1, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesMedianValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesValues = matrix(nrow=nberOfExtractionFiles, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesMeanValue = matrix(nrow=1, ncol=dispersalVelocitySlices)
			weightedBranchDispersalVelocitiesMedianValue = matrix(nrow=1 ,ncol=dispersalVelocitySlices)
			for (i in 1:dispersalVelocitySlices)
				{
					slicedTimes[1,i] = weightedBranchDispersalVelocityList[[1]][i,1]
					for (t in 1:nberOfExtractionFiles)
						{
							weightedBranchDispersalVelocitiesValues[t,i] = weightedBranchDispersalVelocityList[[t]][i,2]
						}
					quantiles = quantile(weightedBranchDispersalVelocitiesValues[,i], probs=c(0.025,0.975), na.rm=T)
					lower_l[1,i] = as.numeric(quantiles[1]); upper_l[1,i] = as.numeric(quantiles[2])
					HPD = HDInterval::hdi(weightedBranchDispersalVelocitiesValues[,i])[1:2]
					lower_l[1,i] = as.numeric(HPD[1]); upper_l[1,i] = as.numeric(HPD[2])
					weightedBranchDispersalVelocitiesMeanValue[1,i] = mean(weightedBranchDispersalVelocitiesValues[,i], na.rm=T)
					weightedBranchDispersalVelocitiesMedianValue[1,i] = median(weightedBranchDispersalVelocitiesValues[,i], na.rm=T)
				}
			yLim = c(0, max(upper_l, na.rm=T))
			treeIDs = paste("distance_tree", treeIDs, sep="")
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = weightedBranchDispersalVelocitiesMeanValue
			tab[,3] = numberOfBranchesMeanValue; colnames(tab) = c("time","velocity","number_of_branches")	
			write.table(tab, file=paste(outputName,"_mean_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab[,1] = slicedTimes; tab[,2] = weightedBranchDispersalVelocitiesMedianValue
			tab[,3] = numberOfBranchesMedianValue; colnames(tab) = c("time","velocity","number_of_branches")	
			write.table(tab, file=paste(outputName,"_median_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")	
			write.table(tab, file=paste(outputName,"_95%HPD_weighted_branch_dispersal_velocity.txt",sep=""), row.names=F, quote=F, sep="\t")

			xLab = "time"; yLab = "weighted branch dispersal velocity"
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_1.pdf",sep=""), width=5, height=5)
			text = "Evolution of weighted branch dispersal velocity"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(weightedBranchDispersalVelocityList[[1]][,1], weightedBranchDispersalVelocityList[[1]][,2], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
			if (nberOfExtractionFiles > 1)
				{
					for (t in 2:nberOfExtractionFiles)
						{
							lines(weightedBranchDispersalVelocityList[[t]][,1],weightedBranchDispersalVelocityList[[t]][,2],lwd=LWD)
						}
				}
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_1.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_weighted_branch_dispersal_velocity_2.pdf",sep=""), width=5, height=5)
			text = "Evolution of weighted branch dispersal velocity"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(slicedTimes, weightedBranchDispersalVelocitiesMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
			xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
			getOption("scipen"); opt = options("scipen"=20)
			polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
			lines(slicedTimes, weightedBranchDispersalVelocitiesMedianValue, lwd=1)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_weighted_branch_dispersal_velocity_2.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
		}
}
