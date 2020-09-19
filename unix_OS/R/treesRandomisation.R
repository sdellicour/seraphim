treesRandomisation <-
function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), randomProcedure=3, nberOfCores=1, showingPlots=F) {

\tnberOfRandomisations = 1
\tregisterDoMC(cores=nberOfCores)
\thull_polygons = list()
\trotation = function(pt1, pt2, angle)
\t\t{
\t\t\ts = sin(angle); c = cos(angle)
\t\t\tx = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
\t\t\tx_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
\t\t\tx_new = x_new+pt1[1]; y_new = y_new+pt1[2]
\t\t\treturn(c(x_new,y_new))
\t\t}
\tnullRaster = envVariables[[1]]
\tnullRaster[!is.na(nullRaster)] = 1
\tnames(nullRaster) = "null_raster"
\tenvVariables0 = envVariables\t
\tnewEnvVariables = list(nullRaster)
\tfor (h in 1:length(envVariables0))
\t\t{
\t\t\tnewEnvVariables[[h+1]] = envVariables0[[h]]
\t\t\tnewEnvVariables[[h+1]][newEnvVariables[[h+1]]<0] = NA
\t\t}
\tenvVariables = newEnvVariables
\tbranchRandomisation3 = FALSE
\tbranchRandomisation2 = FALSE
\tbranchRandomisation1 = FALSE
\tif (randomProcedure == 3) branchRandomisation3 = TRUE
\tif (randomProcedure == 4)
\t\t{
\t\t\tbranchRandomisation2 = TRUE; rotatingEndNodes = TRUE
\t\t}
\tif (randomProcedure == 5)
\t\t{
\t\t\tbranchRandomisation2 = TRUE; rotatingEndNodes = FALSE
\t\t}
\tif (randomProcedure == 6) branchRandomisation1 = TRUE
\tnberOfConnections = matrix(nrow=1, ncol=nberOfExtractionFiles)\t
\ttotalnberOfConnections = 0
\textractionFileName = "TreeExtractions"
\tif (nchar(localTreesDirectory) == 0)
\t\t{
\t\t\tdata = read.csv(paste(extractionFileName,"_1.csv",sep=""), header=T, dec=".")
\t\t}\telse\t{
\t\t\tdata = read.csv(paste(localTreesDirectory,"/",extractionFileName,"_1.csv",sep=""), header=T, dec=".")\t
\t\t}
\tdata = data[with(data, order(startYear,endYear)),]
\tnode1 = list()
\tnode2 = list()
\tstartYear = list()
\tdispersalTime = list()
\ttreeIDs = list()
\tdispersalRate = list()
\tfromCoor = list()
\ttoCoor = list()
\tdatas = list()
\tfor (t in 1:nberOfExtractionFiles)
\t\t{
\t\t\tif (t != 1)
\t\t\t\t{
\t\t\t\t\tif (nchar(localTreesDirectory) == 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfileName = paste(extractionFileName,"_",t,".csv",sep="")
\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\tfileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
\t\t\t\t\t\t}\t
\t\t\t\t\tdata = read.csv(fileName, h = T)
\t\t\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\t\t}
\t\t\tancestralNodeNAonNullRaster = TRUE
\t\t\twhile (ancestralNodeNAonNullRaster == TRUE)
\t\t\t\t{
\t\t\t\t\tancestralNodeNAonNullRaster = FALSE
\t\t\t\t\tancestralBranches = which(!data[,"node1"]%in%data[,"node2"])
\t\t\t\t\tindicesOfBranchesToRemove = c()
\t\t\t\t\tfor (i in 1:length(ancestralBranches))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (is.na(raster::extract(nullRaster, cbind(data[ancestralBranches[i],"startLon"],data[ancestralBranches[i],"startLat"]))))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tancestralNodeNAonNullRaster = TRUE
\t\t\t\t\t\t\t\t\tindicesOfBranchesToRemove = c(indicesOfBranchesToRemove, ancestralBranches[i])
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (length(indicesOfBranchesToRemove) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tdata = data[-indicesOfBranchesToRemove,]
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tnberOfConnections[t] = dim(data)[1]\t\t
\t\t\tnode1[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\tnode1[[t]][] = data[,"node1"]
\t\t\tnode2[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\tnode2[[t]][] = data[,"node2"]
\t\t\tstartYear[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\tstartYear[[t]][] = data[,"startYear"]
\t\t\tdispersalTime[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)
\t\t\tdispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
\t\t\tcolnames(dispersalTime[[t]]) = "dispersalTime"
\t\t\tfromCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\t\tfromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
\t\t\ttoCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\t\ttoCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
\t\t\ttotalnberOfConnections = totalnberOfConnections + nberOfConnections[t] 
\t\t\tif (("treeID"%in%colnames(data)) == TRUE)
\t\t\t\t{
\t\t\t\t\ttreeIDs[[t]] = data[1,"treeID"]
\t\t\t\t}\telse\t{
\t\t\t\t\ttreeIDs[[t]] = "noTreeID"
\t\t\t\t}
\t\t\tdatas[[t]] = data
\t\t}
\thullRasters = list()
\thullRasters[1:length(envVariables)] = envVariables[1:length(envVariables)]
\tpoints = matrix(nrow=(totalnberOfConnections*2),ncol=2)
\ta = 0
\tfor (t in 1:nberOfExtractionFiles)
\t\t{
\t\t\tif (t > 1)
\t\t\t\t{
\t\t\t\t\ta = a + nberOfConnections[t-1]
\t\t\t\t}
\t\t\tfor (i in 1:nberOfConnections[t])
\t\t\t\t{
\t\t\t\t\tindex = (a*2) + ((i-1)*2) + 1
\t\t\t\t\tpoints[index,1] = fromCoor[[t]][i,1]
\t\t\t\t\tpoints[index,2] = fromCoor[[t]][i,2]
\t\t\t\t\tpoints[(index+1),1] = toCoor[[t]][i,1]
\t\t\t\t\tpoints[(index+1),2] = toCoor[[t]][i,2]\t\t\t\t\t\t\t\t
\t\t\t\t}
\t\t}
\tpoints = points[points[,1]>extent(hullRasters[[1]])@xmin,]
\tpoints = points[points[,1]<extent(hullRasters[[1]])@xmax,]
\tpoints = points[points[,2]>extent(hullRasters[[1]])@ymin,]
\tpoints = points[points[,2]<extent(hullRasters[[1]])@ymax,]
\tif (length(hull_polygons) == 0)
\t\t{\t\t
\t\t\thull = chull(points)
\t\t\thull = c(hull,hull[1])
\t\t\t# plot(points); lines(points[hull,])
\t\t\tp = Polygon(points[hull,])
\t\t\tps = Polygons(list(p),1)
\t\t\tsps = SpatialPolygons(list(ps))
\t\t\t# c = mask(envVariables[[2]],sps); plot(c)
\t\t}
\tif (length(hull_polygons) > 0)
\t\t{
\t\t\tsps = hull_polygons
\t\t}
\tpointsRaster = rasterize(points, crop(hullRasters[[1]], sps, snap="out"))
\tpointsRaster[!is.na(pointsRaster[])] = 0\t\t
\tfor (h in 1:length(envVariables))
\t\t{
\t\t\t# plot(mask(simRasters[[h]],sps))
\t\t\thullRasters[[h]] = crop(hullRasters[[h]], sps, snap="out")
\t\t\tbufferRaster = hullRasters[[h]]
\t\t\thullRasters[[h]] = mask(hullRasters[[h]], sps, snap="out")
\t\t\thullRasters[[h]][!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
\t\t\tnames(hullRasters[[h]]) = gsub(".asc","",names(envVariables[[h]]))
\t\t\tnames(hullRasters[[h]]) = gsub(".tif","",names(envVariables[[h]]))
\t\t\tnames(hullRasters[[h]]) = gsub(".gri","",names(envVariables[[h]]))
\t\t}
\tif (nberOfRandomisations > 0)
\t\t{\t
\t\t\tfor (s in 1:nberOfRandomisations)
\t\t\t\t{
\t\t\t\t\tif (branchRandomisation3 == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsimRasters = list(); simRasters = hullRasters
\t\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\t\tbuffer = list()
\t\t\t\t\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t\t\t\t\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\t\t\t\t\tif (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tcounter1 = 0; # print(t)
\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE; counter1 = counter1+1
\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter1 == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcat("Randomising tree ",t,"\\n",sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcat("Randomising tree ",t,", again","\\n",sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][,] = NA; toCoorRand[[t]][,] = NA
\t\t\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# for (i in 1:dim(fromCoor[[t]])[1])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# segments(fromCoor[[t]][i,1], fromCoor[[t]][i,2], toCoor[[t]][i,1], toCoor[[t]][i,2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tancestralIndex = list(); ancestralNodes = list(); counter = 0
\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(node1[[t]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeBoolean = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(node2[[t]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (node1[[t]][i,1] == node2[[t]][j,1])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeBoolean = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeBoolean == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralIndex[[counter]] = i
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodes[[counter]] = node1[[t]][i,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(ancestralIndex))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][ancestralIndex[[i]],1] = fromCoor[[t]][ancestralIndex[[i]],1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][ancestralIndex[[i]],2] = fromCoor[[t]][ancestralIndex[[i]],2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodes = unique(ancestralNodes)
\t\t\t\t\t\t\t\t\t\t\t\t\tstartingNodes = list(); startingNodes = ancestralNodes # startingNodes[[1]] = ancestralNode[1]
\t\t\t\t\t\t\t\t\t\t\t\t\twhile (length(startingNodes) > 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewStartingNodes = list(); c = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(startingNodes))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnodes2 = node2[[t]][which(node1[[t]][,1]==startingNodes[[i]]),1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (length(nodes2) > 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(nodes2))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewStartingNodes[[c]] = nodes2[j]\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tk = which(node2[[t]][,1]==nodes2[j])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt01 = c(fromCoor[[t]][k,1], fromCoor[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt02 = c(toCoor[[t]][k,1], toCoor[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = pt02[1]-pt01[1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = pt02[2]-pt01[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoorRand[[t]][k,1], fromCoorRand[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(NA,NA)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[1] = pt1[1]+xTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[2] = pt1[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter2 = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter2 = counter2+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter2 > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print("counter2 > 100")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeID = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(ancestralNodes))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodes[[h]] == node1[[t]][which(node2[[t]][,1]==nodes2[j]),1])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeID = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeID == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter1 <= 10)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeID == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2[1], pt2[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][k,1] = pt1[1]; fromCoorRand[[t]][k,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][k,1] = pt2[1]; toCoorRand[[t]][k,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoModify = which(node1[[t]][,1]==nodes2[j])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (k in 1:length(toModify))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][toModify[k],1] = pt2[1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][toModify[k],2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tstartingNodes = newStartingNodes\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\ttemp = datas[[t]]
\t\t\t\t\t\t\t\t\t\t\ttemp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
\t\t\t\t\t\t\t\t\t\t\twrite.csv(temp, paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tt
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (branchRandomisation2 == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsimRasters = list(); simRasters = hullRasters
\t\t\t\t\t\t\tcat("Analysis of randomised branch positions ", s, "\\n", sep="")
\t\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\t\tbuffer = list()
\t\t\t\t\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t\t\t\t\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\t\t\t\t\tif (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
\t\t\t\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:nberOfConnections[t])
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif (rotatingEndNodes == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tif (rotatingEndNodes == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt1[1]; toCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt2[1]; fromCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\ttemp = datas[[t]]
\t\t\t\t\t\t\t\t\t\t\ttemp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
\t\t\t\t\t\t\t\t\t\t\twrite.csv(temp, paste0(localTreeDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tt
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}\t\t\t
\t\t\t\t\tif (branchRandomisation1 == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsimRasters = list(); simRasters = hullRasters
\t\t\t\t\t\t\tcat("Analysis of randomised branch positions ", s, "\\n", sep="")
\t\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\t\tbuffer = list()
\t\t\t\t\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t\t\t\t\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\t\t\t\t\tif (!file.exists(paste0(localTreesDirectory,"/TreeRandomisation_",t,".csv")))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\t\t\t\tplot(sps, lwd=0.5, border="black", add=T)
\t\t\t\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:nberOfConnections[t])
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\twhile (pt1NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = runif(1)*(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = runif(1)*(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[1] = pt1[1]+xTranslation; pt1_translated[2] = pt1[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt1_translated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[1] = pt1_translated[1]-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = xTranslation-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt1_translated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[2] = pt1_translated[2]-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = yTranslation-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt1_translated[1],pt1_translated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tinsideAtLeastOneHullPolygon = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(sps))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt1_translated[1], pt1_translated[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tinsideAtLeastOneHullPolygon = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (insideAtLeastOneHullPolygon == FALSE) NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = pt1_translated
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[1] = pt2[1]+xTranslation; pt2[2] = pt2[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_index = NA
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(sps))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt1[1], pt1[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_index = j
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt2_rotated[1], pt2_rotated[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1_translated[1], pt1_translated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\ttemp = datas[[t]]
\t\t\t\t\t\t\t\t\t\t\ttemp[,c("startLon","startLat")] = fromCoorRand[[t]]; temp[,c("endLon","endLat")] = toCoorRand[[t]]
\t\t\t\t\t\t\t\t\t\t\twrite.csv(temp, paste0(localTreeDirectory,"/TreeRandomisation_",t,".csv"), row.names=F, quote=F)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tt
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}\t
\t\t\t\t}
\t\t}\t\t
}
