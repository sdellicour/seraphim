spreadValues <-
function(localTreesDirectory, nberOfExtractionFiles, envVariables, startTime, endTime, timeSlices=200, slidingWindow=1/12, showingPlots=FALSE, outputName=gsub(" ","_",date()), nberOfCores=1, simulations=FALSE) {

\tregisterDoMC(cores=nberOfCores)
\t# for (i in 1:nberOfExtractionFiles)
\t\t# {
\t\t\t# tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
\t\t\t# if (i == 1)
\t\t\t\t# {
\t\t\t\t\t# minStartYear = min(tab[,"startYear"]); maxEndYear = max(tab[,"endYear"])
\t\t\t\t# }\telse\t{
\t\t\t\t\t# if (minStartYear > min(tab[,"startYear"])) minStartYear = min(tab[,"startYear"])
\t\t\t\t\t# if (maxEndYear < max(tab[,"endYear"])) maxEndYear = max(tab[,"endYear"])
\t\t\t\t# }
\t\t# }
\tminStartYear = startTime; maxEndYear = endTime
\ttimeInterval = (maxEndYear-minStartYear)/timeSlices
\tstartEndTimes = matrix(nrow=timeSlices, ncol=3)
\tfor (i in 1:timeSlices)
\t\t{
\t\t\ttime = minStartYear+((i-1)*timeInterval)+(timeInterval/2)
\t\t\tstartTime = time - (slidingWindow/2)
\t\t\tendTime = time + (slidingWindow/2)
\t\t\tstartEndTimes[i,1:3] = cbind(time, startTime, endTime)
\t\t}
\traster_time_intervals = list()
\tenvironmentalValuesList = list(); buffer = list()
\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t# for (t in 1:nberOfExtractionFiles) {
\t\t\tcat(paste0("Analysing tree ",t,"\\n"))
\t\t\tif (simulations == FALSE)
\t\t\t\t{
\t\t\t\t\tfileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
\t\t\t\t}\telse\t{
\t\t\t\t\tfileName = paste(localTreesDirectory,"/TreeSimulations_",t,".csv",sep="")
\t\t\t\t}
\t\t\tdata = read.csv(fileName, h=T)
\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\tnberOfConnections = dim(data)[1]
\t\t\tenvironmentalValues = matrix(nrow=timeSlices, ncol=1+length(envVariables))
\t\t\tdata = data[order(data[,"endYear"]),]
\t\t\tfor (i in 1:timeSlices)
\t\t\t\t{
\t\t\t\t\tenvValues = list()
\t\t\t\t\ttime = startEndTimes[i,1]
\t\t\t\t\tstartTime = startEndTimes[i,2]
\t\t\t\t\tendTime = startEndTimes[i,3]
\t\t\t\t\tfirstExtraction = TRUE
\t\t\t\t\tfor (j in 1:nberOfConnections)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tbranchInTimeInterval = FALSE
\t\t\t\t\t\t\tif ((data[j,"startYear"]<startTime)&(data[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
\t\t\t\t\t\t\tif ((data[j,"startYear"]>startTime)&(data[j,"startYear"]<endTime)) branchInTimeInterval = TRUE
\t\t\t\t\t\t\tif ((data[j,"endYear"]>startTime)&(data[j,"endYear"]<endTime)) branchInTimeInterval = TRUE
\t\t\t\t\t\t\tif (branchInTimeInterval == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (data[j,"startYear"] > startTime) time1 = data[j,"startYear"]
\t\t\t\t\t\t\t\t\tif (data[j,"startYear"] <= startTime) time1 = startTime
\t\t\t\t\t\t\t\t\tif (data[j,"endYear"] < endTime) time2 = data[j,"endYear"]
\t\t\t\t\t\t\t\t\tif (data[j,"endYear"] >= endTime) time2 = endTime
\t\t\t\t\t\t\t\t\tbranchTimeInInterval = time2 - time1
\t\t\t\t\t\t\t\t\ttimeProportion1 = (time1-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
\t\t\t\t\t\t\t\t\tpointLocation1X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion1)
\t\t\t\t\t\t\t\t\tpointLocation1Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion1)
\t\t\t\t\t\t\t\t\tpointLocation1 = cbind(pointLocation1X, pointLocation1Y)
\t\t\t\t\t\t\t\t\ttimeProportion2 = (time2-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
\t\t\t\t\t\t\t\t\tpointLocation2X = data[j,"startLon"]+((data[j,"endLon"]-data[j,"startLon"])*timeProportion2)
\t\t\t\t\t\t\t\t\tpointLocation2Y = data[j,"startLat"]+((data[j,"endLat"]-data[j,"startLat"])*timeProportion2)
\t\t\t\t\t\t\t\t\tpointLocation2 = cbind(pointLocation2X, pointLocation2Y)
\t\t\t\t\t\t\t\t\tbranchDistInInterval = rdist.earth(pointLocation1, pointLocation2, miles=F, R=NULL)
\t\t\t\t\t\t\t\t\tpoints = rbind(cbind(pointLocation1X,pointLocation1Y), cbind(pointLocation2X,pointLocation2Y))
\t\t\t\t\t\t\t\t\tlines = list(); lines[[1]] = Lines(list(Line(points)), 1); spatialLines = SpatialLines(lines)
\t\t\t\t\t\t\t\t\tfor (k in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (dim(envVariables[[k]])[3] == 1)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\textractions = raster::extract(envVariables[[k]], spatialLines); # print(extractions)
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif (i == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttime_intervals = matrix(nrow=length(names(envVariables[[k]])), ncol=2)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (l in 1:length(names(envVariables[[k]])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttime_intervals[l,1] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttime_intervals[l,2] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[3])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\traster_time_intervals[[k]] = time_intervals
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tindex = which((raster_time_intervals[[k]][,1]<time)&(raster_time_intervals[[k]][,2]>time))
\t\t\t\t\t\t\t\t\t\t\t\t\textractions = raster::extract(envVariables[[k]][[index]], spatialLines); # print(extractions)
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tif (firstExtraction == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tenvValues[[k]] = extractions[[1]][!is.na(extractions[[1]])]
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif ((length(envValues)<k)||(is.null(envValues[[k]])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tenvValues[[k]] = extractions[[1]][!is.na(extractions[[1]])]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tenvValues[[k]] = c(envValues[[k]], extractions[[1]][!is.na(extractions[[1]])])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tfirstExtraction = FALSE
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (time > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tenvironmentalValues[i,1] = time
\t\t\t\t\t\t\tif (length(envValues) > 0)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tfor (j in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tenvironmentalValues[i,1+j] = mean(envValues[[j]], rm.na=T)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\t# buffer[[t]] = environmentalValues
\t\t\tenvironmentalValues
\t\t}
\tfor (t in 1:length(buffer))
\t\t{
\t\t\tenvironmentalValuesList[[t]] = buffer[[t]]
\t\t}
\tfor (k in 1:length(envVariables))
\t\t{
\t\t\tif (dim(envVariables[[k]])[3] == 1) envVariableName = names(envVariables[[k]])
\t\t\tif (dim(envVariables[[k]])[3] >= 2) envVariableName = unlist(strsplit(names(envVariables[[k]]),"_"))[1]
\t\t\tslicedTimes = matrix(nrow=1, ncol=timeSlices)
\t\t\tlower_l = matrix(nrow=1, ncol=timeSlices)
\t\t\tupper_l = matrix(nrow=1, ncol=timeSlices)
\t\t\tenvironmentalValues = matrix(nrow=(nberOfExtractionFiles), ncol=timeSlices)
\t\t\tenvironmentalMeanValue = matrix(nrow=1, ncol=timeSlices)
\t\t\tenvironmentalMedianValue = matrix(nrow=1, ncol=timeSlices)
\t\t\tenvironmentalValues = matrix(nrow=(nberOfExtractionFiles), ncol=timeSlices)
\t\t\tenvironmentalMeanValue = matrix(nrow=1, ncol=timeSlices)
\t\t\tfor (i in 1:timeSlices)
\t\t\t\t{
\t\t\t\t\tslicedTimes[1,i] = environmentalValuesList[[1]][i,1]
\t\t\t\t\tfor (t in 1:length(environmentalValuesList))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tenvironmentalValues[t,i] = environmentalValuesList[[t]][i,1+k]
\t\t\t\t\t\t}
\t\t\t\t\tquantiles = quantile(environmentalValues[,i], probs=c(0.025,0.975), na.rm=T)
\t\t\t\t\tlower_l[1,i] = as.numeric(quantiles[1])
\t\t\t\t\tupper_l[1,i] = as.numeric(quantiles[2])
\t\t\t\t\tenvironmentalMeanValue[1,i] = mean(environmentalValues[,i], na.rm=T)
\t\t\t\t\tenvironmentalMedianValue[1,i] = median(environmentalValues[,i], na.rm=T)
\t\t\t\t}
\t\t\txLim = c(minStartYear, maxEndYear); yLim = c(min(upper_l, na.rm=T), max(upper_l, na.rm=T))
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=2)
\t\t\ttab[,1] = slicedTimes; tab[,2] = environmentalMeanValue; colnames(tab) = c("time",envVariableName)\t
\t\t\twrite.csv(tab, file=paste(outputName,"_mean_",envVariableName,".csv",sep=""), row.names=F, quote=F)
\t\t\ttab[,1] = slicedTimes; tab[,2] = environmentalMedianValue; colnames(tab) = c("time",envVariableName)\t
\t\t\twrite.csv(tab, file=paste(outputName,"_median_",envVariableName,".csv",sep=""), row.names=F, quote=F)
\t\t\ttab = matrix(nrow=length(slicedTimes), ncol=3)
\t\t\ttab[,1] = slicedTimes; tab[,2] = lower_l; tab[,3] = upper_l; colnames(tab) = c("time","95%HPD_lower_value","95%HPD_higher_value")\t
\t\t\twrite.csv(tab, file=paste(outputName,"_95%HPD_",envVariableName,".csv",sep=""), row.names=F, quote=F)
\t\t\tLWD = 0.2; xLab = "time"; yLab = envVariableName
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_",envVariableName,"_1.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of environmental values associated with branches"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(environmentalValuesList[[1]][,1], environmentalValuesList[[1]][,1+k], type="l", lwd=0.05, axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\tif (length(environmentalValuesList) > 1)
\t\t\t\t{
\t\t\t\t\tfor (t in 2:length(environmentalValuesList))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tlines(environmentalValuesList[[t]][,1],environmentalValuesList[[t]][,1+k],lwd=LWD)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_",envVariableName,"_1.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t\tif (showingPlots) dev.new(width=5, height=5)
\t\t\tif (showingPlots == FALSE) pdf(paste(outputName,"_",envVariableName,"_2.pdf",sep=""), width=5, height=5)
\t\t\ttext = "Evolution of environmental values associated with branches"
\t\t\tpar(mgp=c(1,0.35,0), oma=c(1,1,1.5,3), mar=c(3.3,3.3,2,0))
\t\t\tplot(slicedTimes, environmentalMedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
\t\t\txx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l,rev(upper_l))
\t\t\tgetOption("scipen"); opt = options("scipen"=20)
\t\t\tpolygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
\t\t\tlines(slicedTimes, environmentalMedianValue, lwd=1)
\t\t\taxis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30")
\t\t\taxis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30")
\t\t\ttitle(xlab=xLab, cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
\t\t\ttitle(ylab=yLab, cex.lab=0.7, mgp=c(1.5,0,0), col.lab="gray30")
\t\t\ttitle(main=text, cex.main=0.6, col.main="gray30")
\t\t\tbox(lwd=LWD, col="gray30")
\t\t\tif (showingPlots) dev.copy2pdf(file=paste(outputName,"_",envVariableName,"_2.pdf",sep=""))
\t\t\tif (showingPlots == FALSE) dev.off()
\t\t}
}
