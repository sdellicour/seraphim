spreadGraphic2 <-
function(localTreesDirectory, nberOfExtractionFiles, prob=0.95, startDatum, precision=1, includeRoot=T) {

\t# nberOfCores=1; registerDoMC(cores=nberOfCores)
\tpercentage = gsub("0\\\\.","",as.character(prob))
\ttimeInterval = precision; nodes = c()
\tfor (i in 1:nberOfExtractionFiles)
\t\t{
\t\t\ttab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
\t\t\tstartingNodeID = which(!tab[,"node1"]%in%tab[,"node2"])[1]
\t\t\tstartingNode = tab[startingNodeID,c("startYear","startLon","startLat")]; colnames(startingNode) = c("time","lon","lat")
\t\t\tendingNodes = tab[,c("endYear","endLon","endLat")]; colnames(endingNodes) = c("time","lon","lat")
\t\t\tif (includeRoot == TRUE) nodes = rbind(nodes, startingNode, endingNodes)
\t\t\tif (includeRoot == FALSE) nodes = rbind(nodes, endingNodes)
\t\t\tif (i == 1) endDatum = max(tab[,"endYear"])
\t\t}
\ttimeSlices = seq(startDatum,endDatum,timeInterval)
\tspreads = list(); c = 0; # buffer = list(); 
\t# buffer = foreach(i = 1:length(timeSlices)) %dopar% {
\tfor (i in 1:length(timeSlices)) {
\t\t\tstartTime = timeSlices[i]-timeInterval/2
\t\t\tendTime = timeSlices[i]+timeInterval/2
\t\t\tif (i == 1) startTime = -9999
\t\t\tselectedNodes = nodes[which((nodes[,"time"]>=startTime)&(nodes[,"time"]<endTime)),]
\t\t\tselectedNodes = unique(selectedNodes[,c("lon","lat")])
\t\t\tif (dim(selectedNodes)[1] > 2)
\t\t\t\t{
\t\t\t\t\tc = c+1; H = Hpi(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"])); print(timeSlices[i])
\t\t\t\t\tkde = kde(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"]), H=H, compute.cont=T, gridsize=c(1000,1000))
\t\t\t\t\tcontourLevel = contourLevels(kde, prob=0.01); polygons = list()
\t\t\t\t\tcontourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
\t\t\t\t\tfor (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
\t\t\t\t\tps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps)); # plot(contourPolygons, add=T)
\t\t\t\t\tspdf = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
\t\t\t\t\tnames(spdf) = round(timeSlices[i],3); spreads[[c]] = spdf; # print(c(i,c,length(spreads)))
\t\t\t\t}
\t\t}
\treturn(spreads)
}
