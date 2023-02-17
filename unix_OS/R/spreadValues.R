spreadValues = function(localTreesDirectory, nberOfExtractionFiles, envVariables, startTime, endTime, timeSlices=200, slidingWindow=1/12, showingPlots=FALSE, outputName=gsub(" ","_",date()), nberOfCores=1, simulations=FALSE) {

	registerDoMC(cores=nberOfCores)
	# for (i in 1:nberOfExtractionFiles)
		# {
			# tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
			# if (i == 1)
				# {
					# minStartYear = min(tab[,"startYear"]); maxEndYear = max(tab[,"endYear"])
				# }	else	{
					# if (minStartYear > min(tab[,"startYear"])) minStartYear = min(tab[,"startYear"])
					# if (maxEndYear < max(tab[,"endYear"])) maxEndYear = max(tab[,"endYear"])
				# }
		# }
	minStartYear = startTime; maxEndYear = endTime
	timeInterval = (maxEndYear-minStartYear)/timeSlices
	startEndTimes = matrix(nrow=timeSlices, ncol=3)
	for (i in 1:timeSlices)
		{
			time = minStartYear+((i-1)*timeInterval)+(timeInterval/2)
			startTime = time - (slidingWindow/2)
			endTime = time + (slidingWindow/2)
			startEndTimes[i,1:3] = cbind(time, startTime, endTime)
		}
	raster_time_intervals = list()
	environmentalValuesList = list(); buffer = list()
	buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
	# for (t in 1:nberOfExtractionFiles) {
			cat(paste0("Analysing tree ",t,"\n"))
			if (simulations == FALSE)
				{
					fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
				}	else	{
					fileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
				}
			data = read.csv(fileName, h=T)
			data = data[with(data, order(endYear, startYear)),]
			nberOfConnections = dim(data)[1]
			environmentalValues = matrix(nrow=timeSlices, ncol=1+length(envVariables))
			data = data[order(data[,"endYear"]),]
			for (i in 1:timeSlices)
				{
					envValues = list()
					time = startEndTimes[i,1]
					startTime = startEndTimes[i,2]
					endTime = startEndTimes[i,3]
					firstExtraction = TRUE
					for (j in 1:nberOfConnections)
						{
							branchInTimeInterval = FALSE
							if ((data[j,"startYear"]<startTime)&(data[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
							if ((data[j,"startYear"]>startTime)&(data[j,"startYear"]<endTime)) branchInTimeInterval = TRUE
							if ((data[j,"endYear"]>startTime)&(data[j,"endYear"]<endTime)) branchInTimeInterval = TRUE
							if (branchInTimeInterval == TRUE)
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
									points = rbind(cbind(pointLocation1X,pointLocation1Y), cbind(pointLocation2X,pointLocation2Y))
									lines = list(); lines[[1]] = Lines(list(Line(points)), 1); spatialLines = SpatialLines(lines)
									for (k in 1:length(envVariables))
										{
											if (dim(envVariables[[k]])[3] == 1)
												{
													extractions = raster::extract(envVariables[[k]], spatialLines); # print(extractions)
												}	else	{
													if (i == 1)
														{
															time_intervals = matrix(nrow=length(names(envVariables[[k]])), ncol=2)
															for (l in 1:length(names(envVariables[[k]])))
																{
																	time_intervals[l,1] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[2])
																	time_intervals[l,2] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[3])
																}
															raster_time_intervals[[k]] = time_intervals
														}
													index = which((raster_time_intervals[[k]][,1]<time)&(raster_time_intervals[[k]][,2]>time))
													extractions = raster::extract(envVariables[[k]][[index]], spatialLines); # print(extractions)
												}
											if (firstExtraction == TRUE)
												{
													envValues[[k]] = extractions[[1]][!is.na(extractions[[1]])]
												}	else	{
													if ((length(envValues)<k)||(is.null(envValues[[k]])))
														{
															envValues[[k]] = extractions[[1]][!is.na(extractions[[1]])]
														}	else	{
															envValues[[k]] = c(envValues[[k]], extractions[[1]][!is.na(extractions[[1]])])
														}
												}
										}
									firstExtraction = FALSE
								}
						}
					if (time > 0)
						{
							environmentalValues[i,1] = time
							if (length(envValues) > 0)
								{
									for (j in 1:length(envVariables))
										{
											environmentalValues[i,1+j] = mean(envValues[[j]], rm.na=T)
										}
								}
						}
				}
			# buffer[[t]] = environmentalValues
			environmentalValues
		}
	for (t in 1:length(buffer))
		{
			environmentalValuesList[[t]] = buffer[[t]]
		}
	for (k in 1:length(envVariables))
		{
			if (dim(envVariables[[k]])[3] == 1) envVariableName = names(envVariables[[k]])
			if (dim(envVariables[[k]])[3] >= 2) envVariableName = unlist(strsplit(names(envVariables[[k]]),"_"))[1]
			slicedTimes = matrix(nrow=1, ncol=timeSlices)
			lower_l = matrix(nrow=1, ncol=timeSlices)
			upper_l = matrix(nrow=1, ncol=timeSlices)
			environmentalValues = matrix(nrow=(nberOfExtractionFiles), ncol=timeSlices)
			environmentalMeanValue = matrix(nrow=1, ncol=timeSlices)
			environmentalMedianValue = matrix(nrow=1, ncol=timeSlices)
			environmentalValues = matrix(nrow=(nberOfExtractionFiles), ncol=timeSlices)
			environmentalMeanValue = matrix(nrow=1, ncol=timeSlices)
			for (i in 1:timeSlices)
				{
					slicedTimes[1,i] = environmentalValuesList[[1]][i,1]
					for (t in 1:length(environmentalValuesList))
						{
							environmentalValues[t,i] = environmentalValuesList[[t]][i,1+k]
						}
					quantiles = quantile(environmentalValues[,i], probs=c(0.025,0.975), na.rm=T)
					lower_l[1,i] = as.numeric(quantiles[1]); upper_l[1,i] = as.numeric(quantiles[2])
					HPD = HDInterval::hdi(environmentalValues[which(!is.na(environmentalValues[,i])),i])[1:2]
					lower_l[1,i] = as.numeric(HPD[1]); upper_l[1,i] = as.numeric(HPD[2])
					environmentalMeanValue[1,i] = mean(environmentalValues[,i], na.rm=T)
					environmentalMedianValue[1,i] = median(environmentalValues[,i], na.rm=T)
				}
			xLim = c(minStartYear, maxEndYear); yLim = c(min(upper_l, na.rm=T), max(upper_l, na.rm=T))
			tab = matrix(nrow=length(slicedTimes), ncol=2)
			tab[,1] = slicedTimes; tab[,2] = environmentalMeanValue; colnames(tab) = c("time",envVariableName)	
			write.csv(tab, file=paste(outputName,"_mean_",envVariableName,".csv",sep=""), row.names=F, quote=F)
			tab[,1] = slicedTimes; tab[,2] = environmentalMedianValue; colnames(tab) = c("time",envVariableName)	
			write.csv(tab, file=paste(outputName,"_median_",envVariableName,".csv",sep=""), row.names=F, quote=F)
			tab = matrix(nrow=length(slicedTimes), ncol=3)
			tab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")	
			write.csv(tab, file=paste(outputName,"_95%HPD_",envVariableName,".csv",sep=""), row.names=F, quote=F)
			LWD = 0.2; xLab = "time"; yLab = envVariableName
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_",envVariableName,"_1.pdf",sep=""), width=5, height=5)
			text = "Evolution of environmental values associated with branches"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(environmentalValuesList[[1]][,1], environmentalValuesList[[1]][,1+k], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
			if (length(environmentalValuesList) > 1)
				{
					for (t in 2:length(environmentalValuesList))
						{
							lines(environmentalValuesList[[t]][,1],environmentalValuesList[[t]][,1+k],lwd=LWD)
						}
				}
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_",envVariableName,"_1.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
			if (showingPlots) dev.new(width=5, height=5)
			if (showingPlots == FALSE) pdf(paste(outputName,"_",envVariableName,"_2.pdf",sep=""), width=5, height=5)
			text = "Evolution of environmental values associated with branches"
			par(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
			plot(slicedTimes, environmentalMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
			xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
			getOption("scipen"); opt = options("scipen"=20)
			polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
			lines(slicedTimes, environmentalMedianValue, lwd=1)
			axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.tick="gray30", col.axis="gray30", col="gray30")
			axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
			title(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
			title(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
			title(main=text, cex.main=0.6, col.main="gray30")
			box(lwd=LWD, col="gray30")
			if (showingPlots) dev.copy2pdf(file=paste(outputName,"_",envVariableName,"_2.pdf",sep=""))
			if (showingPlots == FALSE) dev.off()
		}
}
