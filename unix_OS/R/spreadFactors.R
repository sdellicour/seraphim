spreadFactors <-
function(localTreesDirectory="", nberOfExtractionFiles=1, envVariables=list(), pathModel=1, resistances=list(), avgResistances=list(), fourCells=FALSE, 
\t\t\t nberOfRandomisations=0, randomProcedure=3, outputName="", showingPlots=FALSE, nberOfCores=1, OS="Unix", juliaCSImplementation=FALSE, 
\t\t\t simulations=FALSE, randomisations=FALSE, hull_polygons=list(), onlyTipBranches=FALSE, GLM=FALSE, alternativeQstat=FALSE) {

# variogramModels = list(); hull_polygons = list(); simulations = FALSE; randomisations=FALSE; onlyTipBranches = FALSE; GLM = FALSE; alternativeQstat = FALSE
CA = FALSE
onlyTipBranches = FALSE
impactOnVelocity = TRUE
impactOnDirection = FALSE
if ((pathModel == -1)|(pathModel == 0))
\t{
\t\timpactOnVelocity = FALSE
\t\timpactOnDirection = TRUE
\t}
registerDoMC(cores=nberOfCores)
nberOfCores_CS = 1
plottingHistograms = FALSE
commonalityAnalysis = FALSE
if (CA == TRUE) commonalityAnalysis = TRUE
all = FALSE # boolean variable specifying if LR coefficients have be tested with randomisation on all rasters (TRUE) or with one environmental 
\t# raster at a time (FALSE). In the case of torus translations, if TRUE, all the rasters will be uniformly translated at each iteration
thetaValue = 1 # numeric value [0 < theta < 20], the degree from which the path randomly deviates from the shortest path
envVariableToRandomise = -1 # index of the environmental data to randomise. By default (= -1), all the data will be randomised
dispersalTimeBoolean = TRUE
date1 = base::date()
preLogTransformation = function(x, m)
\t{
\t\tx = 1+(x-m)
\t}
logTransformation = function(x, m)
\t{
\t\tx = log10(x)
\t}
zTransformation = function(x)
\t{ 
\t\tx = (x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T))
\t}
featureScaling = function(x)
\t{ 
\t\tx = (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
\t}\t
rotation = function(pt1, pt2, angle)
\t{
\t\ts = sin(angle); c = cos(angle)
\t\tx = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
\t\tx_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
\t\tx_new = x_new+pt1[1]; y_new = y_new+pt1[2]
\t\treturn(c(x_new,y_new))
\t}
nullRaster = envVariables[[1]]
nullRaster[!is.na(nullRaster)] = 1
names(nullRaster) = "null_raster"
envVariables0 = envVariables\t
newEnvVariables = list(nullRaster)
newResistances = c(TRUE)
newAvgResistances = c(TRUE)
for (h in 1:length(envVariables0))
\t{
\t\tnewEnvVariables[[h+1]] = envVariables0[[h]]
\t\tnewEnvVariables[[h+1]][newEnvVariables[[h+1]]<0] = NA
\t\tif(resistances[[h]] == TRUE)
\t\t\t{
\t\t\t\tfric = "_R"
\t\t\t}\telse\t\t{
\t\t\t\tfric = "_C"
\t\t\t}
\t\tnames(newEnvVariables[[h+1]]) = paste(names(newEnvVariables[[h+1]]), fric, sep="")
\t\tif (length(resistances) > 0)
\t\t\t{
\t\t\t\tnewResistances = c(newResistances, resistances[[h]])
\t\t\t\tnewAvgResistances = c(newAvgResistances, avgResistances[[h]])
\t\t\t}
\t}
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
\t{
\t\tbranchRandomisation2 = TRUE; rotatingEndNodes = TRUE
\t}
if (randomProcedure == 5)
\t{
\t\tbranchRandomisation2 = TRUE; rotatingEndNodes = FALSE
\t}
if (randomProcedure == 6) branchRandomisation1 = TRUE
if (randomProcedure == 7) distPermutations = TRUE
nberOfConnections = matrix(nrow=1, ncol=nberOfExtractionFiles)\t
totalnberOfConnections = 0
if ((simulations == FALSE)&(randomisations == FALSE))
\t{
\t\textractionFileName = "TreeExtractions"
\t}
if ((simulations == TRUE)&(randomisations == FALSE))
\t{
\t\textractionFileName = "TreeSimulations"
\t}
if ((simulations == FALSE)&(randomisations == TRUE))
\t{
\t\textractionFileName = "TreeRandomisation"
\t}
if (nchar(localTreesDirectory) == 0)
\t{
\t\tdata = read.csv(paste(extractionFileName,"_1.csv",sep=""), header=T, dec=".")
\t}\telse\t\t{
\t\tdata = read.csv(paste(localTreesDirectory,"/",extractionFileName,"_1.csv",sep=""), header=T, dec=".")\t
\t}
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
\t{
\t\tdistances = list()
\t}
if (impactOnDirection == TRUE)
\t{
\t\tmeanEnvValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\tmeanDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t}
for (t in 1:nberOfExtractionFiles)
\t{
\t\tif (t != 1)
\t\t\t{
\t\t\t\tif (nchar(localTreesDirectory) == 0)
\t\t\t\t\t{
\t\t\t\t\t\tfileName = paste(extractionFileName,"_",t,".csv",sep="")
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tfileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
\t\t\t\t\t}\t
\t\t\t\tdata = read.csv(fileName, h = T)
\t\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\t}
\t\tancestralNodeNAonNullRaster = TRUE
\t\twhile (ancestralNodeNAonNullRaster == TRUE)
\t\t\t{
\t\t\t\tancestralNodeNAonNullRaster = FALSE
\t\t\t\tancestralBranches = which(!data[,"node1"]%in%data[,"node2"])
\t\t\t\tindicesOfBranchesToRemove = c()
\t\t\t\tfor (i in 1:length(ancestralBranches))
\t\t\t\t\t{
\t\t\t\t\t\tif (is.na(raster::extract(nullRaster, cbind(data[ancestralBranches[i],"startLon"],data[ancestralBranches[i],"startLat"]))))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tancestralNodeNAonNullRaster = TRUE
\t\t\t\t\t\t\t\tindicesOfBranchesToRemove = c(indicesOfBranchesToRemove, ancestralBranches[i])
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (length(indicesOfBranchesToRemove) > 0)
\t\t\t\t\t{
\t\t\t\t\t\tdata = data[-indicesOfBranchesToRemove,]
\t\t\t\t\t}
\t\t\t}
\t\tif (onlyTipBranches == TRUE)
\t\t\t{
\t\t\t\tindices = which(!data[,"node2"]%in%data[,"node1"])
\t\t\t\tdata = data[indices,]
\t\t\t}
\t\tnberOfConnections[t] = dim(data)[1]\t\t
\t\tnode1[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\tnode1[[t]][] = data[,"node1"]
\t\tnode2[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\tnode2[[t]][] = data[,"node2"]
\t\tstartYear[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\tstartYear[[t]][] = data[,"startYear"]
\t\tdispersalTime[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)
\t\tdispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
\t\tcolnames(dispersalTime[[t]]) = "dispersalTime"
\t\tfromCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\tfromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
\t\ttoCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\ttoCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
\t\ttotalnberOfConnections = totalnberOfConnections + nberOfConnections[t] 
\t\tif (impactOnVelocity == TRUE)
\t\t\t{
\t\t\t\tdistances[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
\t\t\t}
\t\tif (("treeID"%in%colnames(data)) == TRUE)
\t\t\t{
\t\t\t\ttreeIDs[[t]] = data[1,"treeID"]
\t\t\t}\telse\t{
\t\t\t\ttreeIDs[[t]] = "noTreeID"
\t\t\t}\t
\t}
hullRasters = list()
hullRasters[1:length(envVariables)] = envVariables[1:length(envVariables)]
points = matrix(nrow=(totalnberOfConnections*2),ncol=2)
a = 0
for (t in 1:nberOfExtractionFiles)
\t{
\t\tif (t > 1)
\t\t\t{
\t\t\t\ta = a + nberOfConnections[t-1]
\t\t\t}
\t\tfor (i in 1:nberOfConnections[t])
\t\t\t{
\t\t\t\tindex = (a*2) + ((i-1)*2) + 1
\t\t\t\tpoints[index,1] = fromCoor[[t]][i,1]
\t\t\t\tpoints[index,2] = fromCoor[[t]][i,2]
\t\t\t\tpoints[(index+1),1] = toCoor[[t]][i,1]
\t\t\t\tpoints[(index+1),2] = toCoor[[t]][i,2]\t\t\t\t\t\t\t\t
\t\t\t}
\t}
points = points[points[,1]>extent(hullRasters[[1]])@xmin,]
points = points[points[,1]<extent(hullRasters[[1]])@xmax,]
points = points[points[,2]>extent(hullRasters[[1]])@ymin,]
points = points[points[,2]<extent(hullRasters[[1]])@ymax,]
if (length(hull_polygons) == 0)
\t{\t\t
\t\thull = chull(points)
\t\thull = c(hull,hull[1])
\t\t# plot(points); lines(points[hull,])
\t\tp = Polygon(points[hull,])
\t\tps = Polygons(list(p),1)
\t\tsps = SpatialPolygons(list(ps))
\t\t# c = mask(envVariables[[2]],sps); plot(c)
\t}
if (length(hull_polygons) > 0)
\t{
\t\tsps = hull_polygons
\t}
pointsRaster = rasterize(points, crop(hullRasters[[1]], sps, snap="out"))
pointsRaster[!is.na(pointsRaster[])] = 0\t\t
for (h in 1:length(envVariables))
\t{
\t\t# plot(mask(simRasters[[h]],sps))
\t\thullRasters[[h]] = crop(hullRasters[[h]], sps, snap="out")
\t\tbufferRaster = hullRasters[[h]]
\t\thullRasters[[h]] = raster::mask(hullRasters[[h]], sps, snap="out")
\t\thullRasters[[h]][!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
\t\tnames(hullRasters[[h]]) = gsub(".asc","",names(envVariables[[h]]))
\t\tnames(hullRasters[[h]]) = gsub(".tif","",names(envVariables[[h]]))
\t\tnames(hullRasters[[h]]) = gsub(".gri","",names(envVariables[[h]]))
\t}
if (randomWalkDistance == TRUE)
\t{
\t\textensions = rep("", length(envVariables))
\t\tif ("CS_rasters"%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), "CS_rasters"))
\t\tfor (h in 1:length(envVariables))
\t\t\t{
\t\t\t\textensions[h] = ".asc"
\t\t\t\tif (round(res(hullRasters[[h]])[1],10) != round(res(hullRasters[[h]])[2],10)) extensions[h] = ".tif"
\t\t\t\tname = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
\t\t\t\twriteRaster(hullRasters[[h]], name, overwrite=T)
\t\t\t}
\t}

## 1. Estimation of the different environmental distances or differences associated with each branch

if (distanceMatrix == F)
\t{
\t\tfor (h in 1:length(envVariables))
\t\t\t{
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tif (straightLineDistance == TRUE) cat("Computing environmental distances (straight-line path model) for ", names(envVariables[[h]])[1], "\\n", sep="")
\t\t\t\t\t\tif (leastCostDistance == TRUE) cat("Computing environmental distances (least-cost path model) for ", names(envVariables[[h]])[1], "\\n", sep="")
\t\t\t\t\t\tif (randomWalkDistance == TRUE) cat("Computing environmental distances (Circuitscape path model) for ", names(envVariables[[h]])[1], "\\n", sep="")\t
\t\t\t\t\t\tif (fourCells == TRUE) directions = 4
\t\t\t\t\t\tif (fourCells == FALSE) directions = 8
\t\t\t\t\t\tif (resistances[h] == FALSE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (randomWalkDistance == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttrEnvVariable = transition(hullRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (showingPlots == FALSE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tbuffer = list()
\t\t\t\t\t\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t\t\t\t\t\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\t\t\t\t\t\tmat = matrix(nrow=nberOfConnections[t], ncol=1)
\t\t\t\t\t\t\t\t\t\tif (straightLineDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tlinesList = list()
\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoor[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints = rbind(fromCoor[[t]][i,], toCoor[[t]][i,])
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlinesList[[i]] = Lines(list(Line(points)),i)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tlines = SpatialLines(linesList)
\t\t\t\t\t\t\t\t\t\t\t\textractions = raster::extract(hullRasters[[h]], lines)
\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoor[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[i] = sum(extractions[[i]], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmat[] = diag(costDistance(trEnvVariableCorr, fromCoor[[t]], toCoor[[t]]))\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (randomWalkDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tenvVariableName = paste("CS_rasters/",names(hullRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
\t\t\t\t\t\t\t\t\t\t\t\tbranchesNotNA = which(!((is.na(raster::extract(hullRasters[[h]], fromCoor[[t]][])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  |(is.na(raster::extract(hullRasters[[h]], toCoor[[t]][])))))
\t\t\t\t\t\t\t\t\t\t\t\tfromCoor_temp = fromCoor[[t]][branchesNotNA,]; toCoor_temp = toCoor[[t]][branchesNotNA,]
\t\t\t\t\t\t\t\t\t\t\t\tif (juliaCSImplementation == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[branchesNotNA,1] = circuitScape1(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t# if (-777%in%mat[])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# mat[branchesNotNA,1] = circuitScape1(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t # fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[branchesNotNA,1] = circuitScape2(hullRasters[[h]],envVariableName,resistances[[h]],avgResistances[[h]],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t fourCells,fromCoor_temp,toCoor_temp,OS,outputName,t,nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoor[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tspatialPoints = SpatialPoints(cbind(c(fromCoor[[t]][i,1], toCoor[[t]][i,1]), c(fromCoor[[t]][i,2], toCoor[[t]][i,2])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[i] = commuteDistance(trEnvVariableCorr, spatialPoints)
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmat[] = diag(rSPDistance(trEnvVariableCorr, fromCoor[[t]], toCoor[[t]], theta=thetaValue, totalNet="total", method=1))
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\tcolnames(mat) = names(envVariables[[h]])
\t\t\t\t\t\t\t\t\t\t# buffer[[t]] = mat
\t\t\t\t\t\t\t\t\t\tmat
\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\tfor (t in 1:length(buffer))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tbuffer[[t]][!is.finite(buffer[[t]][])] = NA # least-cost case
\t\t\t\t\t\t\t\t\t\tbuffer[[t]][buffer[[t]][]==-1] = NA # CircuitScape (RW) case
\t\t\t\t\t\t\t\t\t\tbuffer[[t]][buffer[[t]][]==-777] = NA # NA value in Circuitscape
\t\t\t\t\t\t\t\t\t\tdistances[[t]][,h] = buffer[[t]][]
\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tplotRaster(envVariables[[h]], addLegend=T)
\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\ttext1 = paste(names(envVariables[[h]])[1],", sampled tree ",t,sep="")
\t\t\t\t\t\t\t\t\t\tif (straightLineDistance == TRUE) text2 = "Computing straight-line distances..."
\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE) text2 = "Computing least-cost distances..."
\t\t\t\t\t\t\t\t\t\tif ((randomWalkDistance == TRUE) | (commuteDistance == TRUE)) text2 = "Computing Circuitscape distances..."
\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE) text2 = "Computing randomised shortest path distances..."
\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t\tmtext(text2, col="red", cex=0.7, line=-1)
\t\t\t\t\t\t\t\t\t\tfor (j in 1:nberOfConnections[t])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (j == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(fromCoor[[t]][j,1], fromCoor[[t]][j,2], pch=16, col="black", cex=0.5)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tsegments(fromCoor[[t]][j,1], fromCoor[[t]][j,2], toCoor[[t]][j,1], toCoor[[t]][j,2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\tpoints(toCoor[[t]][j,1], toCoor[[t]][j,2], pch=16, col="black", cex=0.5)
\t\t\t\t\t\t\t\t\t\t\t\tif (straightLineDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlinesList = list()
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints = rbind(fromCoor[[t]][j,], toCoor[[t]][j,])
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlinesList[[1]] = Lines(list(Line(points)),j)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlines = SpatialLines(linesList)
\t\t\t\t\t\t\t\t\t\t\t\t\t\textractions = raster::extract(envVariables[[h]], lines)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = sum(extractions[[1]], na.rm=T)\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorJ = cbind(fromCoor[[t]][j,1], fromCoor[[t]][j,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorJ = cbind(toCoor[[t]][j,1], toCoor[[t]][j,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = costDistance(trEnvVariableCorr, fromCoorJ, toCoorJ)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (randomWalkDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tenvVariableName = paste("CS_rasters/",names(envVariables[[h]]),"_",outputName,"_cs",extensions[h],sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromC = matrix(nrow=1, ncol=2); toC = matrix(nrow=1, ncol=2)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromC[] = fromCoor[[t]][j,]; toC[] = toCoor[[t]][j,]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (juliaCSImplementation == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = circuitScape1(envVariables[[h]], envVariableName, resistances[[h]], avgResistances[[h]], 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   \t\t\tfourCells, fromC, toC, OS, outputName,t, nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = circuitScape2(envVariables[[h]], envVariableName, resistances[[h]], avgResistances[[h]], 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   \t\t\tfourCells, fromC, toC, OS, outputName,t, nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tspatialPoints = SpatialPoints(cbind(c(fromCoor[[t]][j,1], toCoor[[t]][j,1]), c(fromCoor[[t]][j,2], toCoor[[t]][j,2])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = commuteDistance(trEnvVariableCorr, spatialPoints)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = rSPDistance(trEnvVariableCorr, fromCoor[[t]][j,], toCoor[[t]][j,],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  theta=thetaValue, totalNet="total", method=1)
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(as.numeric(distances[[t]][j,h])) == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (!is.finite(distances[[t]][j,h]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistances[[t]][j,h] = NA
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (distances[[t]][j,h] == -1) distances[[t]][j,h] = NA\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tdev.off()
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (impactOnDirection == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tenvValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
\t\t\t\t\t\t\t\tfor (i in 1:length(ancestralNodes))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
\t\t\t\t\t\t\t\t\t\tenvValues = envValues + raster::extract(envVariables[[h]], cbind(fromCoor[[t]][ancestralBranch,1],fromCoor[[t]][ancestralBranch,2]))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tmeanEnvValues[t,h] = envValues + mean(raster::extract(envVariables[[h]], toCoor[[t]]), na.rm=T)
\t\t\t\t\t\t\t\tdiffs = raster::extract(envVariables[[h]], fromCoor[[t]])-raster::extract(envVariables[[h]], toCoor[[t]])
\t\t\t\t\t\t\t\trateOfPositiveDifferences[t,h] = sum(diffs[!is.na(diffs)] > 0)/length(diffs[!is.na(diffs)])
\t\t\t\t\t\t\t\tmeanDifferences[t,h] = mean(diffs, na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}\t\t\t\t\t
\t}

envVariableNames = names(envVariables[[1]])\t
for (h in 2:length(envVariables))
\t{
\t\tenvVariableNames = cbind(envVariableNames, names(envVariables[[h]]))
\t}
for (t in 1:nberOfExtractionFiles)
\t{
\t\tif (impactOnVelocity == TRUE) colnames(distances[[t]]) = envVariableNames
\t}

# To generate scatterplots:
if ((showingPlots == TRUE) & (impactOnVelocity == TRUE))
\t{
\t\t# for (h in 1:length(envVariables))
\t\t\t# {
\t\t\t\t# fileName = paste(outputName,"_scatterplot_",names(envVariables[[h]]),".pdf",sep="")
\t\t\t\t# dists = distances[[1]][,h]
\t\t\t\t# times = dispersalTime[[1]]
\t\t\t\t# if (nberOfExtractionFiles > 1)
\t\t\t\t\t# {
\t\t\t\t\t\t# for (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t# dists = c(dists, distances[[t]][,h])
\t\t\t\t\t\t\t\t# times = c(times, dispersalTime[[t]])
\t\t\t\t\t\t\t# }
\t\t\t\t\t# }
\t\t\t\t# dev.new(width=5, height=5); par(mar=c(3.5,3.5,4,2), mgp=c(1.6,0.5,0))\t
\t\t\t\t# formula = paste("times ~ dists", sep="")
\t\t\t\t# form = as.formula(formula); LM = lm(form)
\t\t\t\t# title = paste("Linear regression R2 = ",round(summary(LM)$r.squared,3),sep="")
\t\t\t\t# plot(dists, times, xlab="environmental distances", ylab="dispersal times", main=title, cex.main=0.7, cex.axis=0.6, lwd=0.8, pch=19, cex=0.3, cex.lab=0.7)
\t\t\t\t# abline(LM, lwd=1.5, col="red")
\t\t\t\t# dev.copy2pdf(file=fileName)
\t\t\t\t# dev.off()
\t\t\t# }
\t}
if ((nberOfExtractionFiles == 1) & (impactOnVelocity == TRUE))
\t{
\t\tfileName = paste(outputName, "_env_distances.txt", sep="")
\t\tmat = cbind(dispersalTime[[1]][], distances[[1]][,1:length(envVariables)])
\t\tcolumnNames = "dispersal_times"
\t\tfor (h in 1:length(envVariables))
\t\t\t{
\t\t\t\tcolumnNames = cbind(columnNames, names(envVariables[[h]]))
\t\t\t}
\t\tcolnames(mat) = columnNames
\t\twrite.table(mat, file=fileName, row.names=F, quote=F, sep="\\t")
\t}
if ((file.exists(outputName)) & (impactOnVelocity == TRUE))
\t{\t
\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t{
\t\t\t\tfileName = paste(outputName, "/", outputName, "_tree", t, "_env_distances.txt", sep="")
\t\t\t\tmat = cbind(dispersalTime[[t]][], distances[[t]][,1:length(envVariables)])
\t\t\t\tcolumnNames = cbind("dispersal_times")
\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tcolumnNames = cbind(columnNames, names(envVariables[[h]]))
\t\t\t\t\t}
\t\t\t\tcolnames(mat) = columnNames
\t\t\t\twrite.table(mat, file=fileName, row.names=F, quote=F, sep="\\t")
\t\t\t}
\t}

## 2. Linear regressions of dispersal time vs. environmental distances

if (impactOnVelocity == TRUE)
\t{
\t\t# realUniLRresiduals = list()
\t\trealUniLRcoefficients1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniLRcoefficients2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniLRRsquarePValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\trealUniDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\tif (GLM == TRUE)
\t\t\t{
\t\t\t\trealMultiGLMcoefficients = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trealMultiGLMresiduals = matrix(nrow=length(dispersalTime[[1]]), ncol=nberOfExtractionFiles)
\t\t\t\trealMultiGLMcoefficientPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trealMultiGLMCAuniqueContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trealMultiGLMCAcommonContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t}
\t\tcolNames = c(); for (t in 1:nberOfExtractionFiles) colNames = c(colNames, paste("tree_",treeIDs[[t]],sep=""))
\t\tif (GLM == TRUE) colnames(realMultiGLMresiduals) = colNames
\t\tfor (h in 1:length(envVariables))
\t\t\t{
\t\t\t\tnRowsMax = length(dispersalTime[[1]])
\t\t\t\tif (nberOfExtractionFiles > 1)
\t\t\t\t\t{
\t\t\t\t\t\tfor (t in 2:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (nRowsMax < length(dispersalTime[[t]])) nRowsMax = length(dispersalTime[[t]])
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\t# residuals = matrix(nrow=nRowsMax, ncol=nberOfExtractionFiles)
\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t{
\t\t\t\t\t\tdistVariables = paste("dispersalTime[[",t,"]]"," ~ distances[[",t,"]][,",h,"]",sep="")
\t\t\t\t\t\tform = as.formula(distVariables)
\t\t\t\t\t\tLM = lm(form)
\t\t\t\t\t\trealUniLRcoefficients1[t,h] = summary(LM)$coefficients[2,"Estimate"]
\t\t\t\t\t\trealUniLRRsquares1[t,h] = summary(LM)$r.squared
\t\t\t\t\t\tf = summary(LM)$fstatistic
\t\t\t\t\t\tif (is.numeric(f))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tp = pf(f[1],f[2],f[3],lower.tail=F)
\t\t\t\t\t\t\t\tattributes(p) = NULL; realUniLRRsquarePValues[t,h] = p
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\trealUniLRRsquarePValues[t,h] = NA 
\t\t\t\t\t\t\t}
\t\t\t\t\t\tnonNaNdispersalTimes = dispersalTime[[t]][!is.na(distances[[t]][,1])]
\t\t\t\t\t\tnonNaNdistances1 = distances[[t]][,1][!is.na(distances[[t]][,1])]
\t\t\t\t\t\tnonNaNdistances2 = distances[[t]][,h][!is.na(distances[[t]][,1])]
\t\t\t\t\t\tdistVariables = paste("nonNaNdistances2 ~ nonNaNdistances1",sep="")
\t\t\t\t\t\tLM1 = lm(as.formula(distVariables)); residuals_LM1 = LM1$residuals; indices = which(!is.na(residuals_LM1))
\t\t\t\t\t\tdistVariables = paste("nonNaNdispersalTimes[indices] ~ nonNaNdistances2[indices] + residuals_LM1",sep="")
\t\t\t\t\t\tLM2 = lm(as.formula(distVariables)); realUniLRRsquares2[t,h] = summary(LM2)$r.squared
\t\t\t\t\t\trealUniLRcoefficients2[t,h] = summary(LM2)$coefficients[2,"Estimate"]
\t\t\t\t\t\t# if (length(residuals[,t]) == length(studres(LM)))
\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t# residuals[,t] = studres(LM)
\t\t\t\t\t\t\t# }\telse\t\t{
\t\t\t\t\t\t\t\t# for (i in 1:length(studres(LM)))
\t\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t\t# residuals[i,t] = studres(LM)[i]
\t\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t# }
\t\t\t\t\t}
\t\t\t\t# realUniLRresiduals[[h]] = residuals
\t\t\t}
\t\tif (GLM == TRUE)
\t\t\t{
\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t{
\t\t\t\t\t\tdistVariables = ""
\t\t\t\t\t\tmatMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
\t\t\t\t\t\tmatMultiGLM[,1] = dispersalTime[[t]]
\t\t\t\t\t\tmatMultiGLMNames = c("dispersalTime")
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmatMultiGLM[,h+1] = distances[[t]][,h]
\t\t\t\t\t\t\t\tmatMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
\t\t\t\t\t\t\t\tif (h == 1)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tdistVariables = paste("dispersalTime ~ ", names(envVariables[[h]]), sep="")
\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\tdistVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\tm = min(unlist(lapply(as.matrix(matMultiGLM),min)), na.rm=T)\t
\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, preLogTransformation, m)
\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, logTransformation)
\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, zTransformation)
\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\tform = as.formula(distVariables)
\t\t\t\t\t\tmultiGLM = stats::glm(form, data=matMultiGLM)
\t\t\t\t\t\tif (length(realMultiGLMresiduals[,t]) == length(studres(multiGLM)))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\trealMultiGLMresiduals[,t] = studres(multiGLM)
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tfor (i in 1:length(studres(multiGLM)))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\trealMultiGLMresiduals[i,t] = studres(multiGLM)[i]
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tCA = calc.yhat(multiGLM, prec=5)\t
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tnames(multiGLM$coefficients)[1+h] = names(envVariables[[h]])
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\tname = names(envVariables[[h]]); nameInside = FALSE
\t\t\t\t\t\t\t\tfor (i in 1:length(summary(multiGLM)[]$coefficients[,"Estimate"]))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (name == rownames(summary(multiGLM)[]$coefficients)[i])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t nameInside = TRUE
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (nameInside == TRUE)
\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\trealMultiGLMcoefficients[t,h] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
\t\t\t\t\t\t\t\t\t\trealMultiGLMcoefficientPValues[t,h] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\trealMultiGLMCAuniqueContributions[t,h] = CA$PredictorMetrics[h,"Unique"]
\t\t\t\t\t\t\t\t\t\t\t\trealMultiGLMCAcommonContributions[t,h] = CA$PredictorMetrics[h,"Common"]\t
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}
\t\t\t\t
\t\tuniLRmedianRsquares = matrix(nrow=1, ncol=length(envVariables))
\t\tuniLRmedianRsquarePValues = matrix(nrow=1, ncol=length(envVariables))
\t\tuniLRmedianDeltaRsquares = matrix(nrow=1, ncol=length(envVariables))
\t\tif (GLM == TRUE)
\t\t\t{
\t\t\t\tmultiGLMmedianCoefficients = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tmultiGLMmedianCoefficientPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tmultiGLMmedianCAuniqueContributions = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tmultiGLMmedianCAcommonContributions = matrix(nrow=1, ncol=length(envVariables))
\t\t\t}
\t\tfor (h in 1:length(envVariables))
\t\t\t{
\t\t\t\tuniLRmedianRsquares[1,h] = median(realUniLRRsquares1[,h], na.rm=T)
\t\t\t\tuniLRmedianRsquarePValues[1,h] = median(realUniLRRsquarePValues[,h], na.rm=T)
\t\t\t\tuniLRmedianDeltaRsquares[1,h] = median(realUniLRRsquares1[,h]-realUniLRRsquares1[,1], na.rm=T)
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMmedianCoefficients[1,h] = median(realMultiGLMcoefficients[,h], na.rm=T)
\t\t\t\t\t\tmultiGLMmedianCoefficientPValues[1,h] = median(realMultiGLMcoefficientPValues[,h], na.rm=T)
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMmedianCAuniqueContributions[1,h] = median(realMultiGLMCAuniqueContributions[,h], na.rm=T)
\t\t\t\t\t\t\t\tmultiGLMmedianCAcommonContributions[1,h] = median(realMultiGLMCAcommonContributions[,h], na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (h != 1)
\t\t\t\t\t{
\t\t\t\t\t\trealUniDeltaRsquares1[,h] = (realUniLRRsquares1[,h]-realUniLRRsquares1[,1])
\t\t\t\t\t\trealUniDeltaRsquares2[,h] = (realUniLRRsquares2[,h]-realUniLRRsquares1[,1])
\t\t\t\t\t}
\t\t\t}
\t\t\t
\t\tif ((plottingHistograms == TRUE)&(nberOfExtractionFiles > 1))
\t\t\t{
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfileName1 = paste(outputName, "_LR-GLM-CA_results.pdf", sep="")
\t\t\t\t\t\t\t\tpdf(fileName1, width=(4*(length(envVariables))), height=(5*4)); par(mfrow=c(6,(length(envVariables))))
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tfileName1 = paste(outputName, "_LR-GLM_results.pdf", sep="")
\t\t\t\t\t\t\t\tpdf(fileName1, width=(4*(length(envVariables))), height=(3*4)); par(mfrow=c(4,(length(envVariables))))
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tfileName1 = paste(outputName, "_linear_regression_results.pdf", sep="")
\t\t\t\t\t\tpdf(fileName1, width=(4*(length(envVariables))), height=(2.25*4)); par(mfrow=c(3,(length(envVariables))))
\t\t\t\t\t}
\t\t\t\tbreakList_1 = (0:50)/50; breakList_2 = (0:50)/100
\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\txMin = min(realUniLRcoefficients1[,h]); xMax = max(realUniLRcoefficients1[,h]) # breaks=seq(min,max,by=(max-min)/50)
\t\t\t\t\t\thist(realUniLRcoefficients1[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50), 
\t\t\t\t\t\t\t xlab="Univariate LR coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t}
\t\t\t\txMin = min(realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,1], na.rm=T)\t
\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tif (xMin > min(realUniLRRsquares1[,h], na.rm=T))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMin = min(realUniLRRsquares1[,h], na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (xMax < max(realUniLRRsquares1[,h], na.rm=T))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMax = max(realUniLRRsquares1[,h], na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tmin = min(realUniLRRsquares1[,h]); max = max(realUniLRRsquares1[,h]) # breaks=seq(min,max,by=(max-min)/50)
\t\t\t\t\t\thist(realUniLRRsquares1[,h], freq=T, xlim=c(0,xMax), breaks=seq(0,xMax,by=(xMax-0)/50), xlab="Univariate LR R2's", 
\t\t\t\t\t\t\t main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t}
\t\t\t\tplot.new()
\t\t\t\txMin = min(realUniLRRsquares1[,2]-realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,2]-realUniLRRsquares1[,1], na.rm=T)
\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tif (xMin > min(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMin = min(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (xMax < max(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMax = max(realUniLRRsquares1[,h]-realUniLRRsquares1[,1],na.rm=T)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\thist(realUniDeltaRsquares1[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50), xlab="Univariate LR delta R2 (Q)", 
\t\t\t\t\t\t\t main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t}
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\txMin = min(realMultiGLMcoefficients[,1], na.rm=T); xMax = max(realMultiGLMcoefficients[,1], na.rm=T)
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (xMin > min(realMultiGLMcoefficients[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\txMin = min(realMultiGLMcoefficients[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (xMax < max(realMultiGLMcoefficients[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\txMax = max(realMultiGLMcoefficients[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (is.na(mean(realMultiGLMcoefficients[,h])))
\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\thist(realMultiGLMcoefficients[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}\t
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMin = min(c(realMultiGLMCAuniqueContributions[,1],realMultiGLMCAcommonContributions[,1]), na.rm=T)
\t\t\t\t\t\t\t\txMax = max(c(realMultiGLMCAuniqueContributions[,1],realMultiGLMCAcommonContributions[,1]), na.rm=T)
\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (xMin > min(realMultiGLMCAuniqueContributions[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMin = min(realMultiGLMCAuniqueContributions[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (xMin > min(realMultiGLMCAcommonContributions[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMin = min(realMultiGLMCAcommonContributions[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\tif (xMax < max(realMultiGLMCAuniqueContributions[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMax = max(realMultiGLMCAuniqueContributions[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (xMax < max(realMultiGLMCAcommonContributions[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMax = max(realMultiGLMCAcommonContributions[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (is.na(mean(realMultiGLMcoefficients[,h])))
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\thist(realMultiGLMCAuniqueContributions[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
\t\t\t\t\t\t\t\t\t\t\t\txlab="Multivariate GLM-CA unique contributions", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (is.na(mean(realMultiGLMcoefficients[,h])))
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\thist(realMultiGLMCAcommonContributions[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
\t\t\t\t\t\t\t\t\t\t\t\txlab="Multivariate GLM-CA common contributions", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tdev.off()
\t\t\t}
\t\tif (nberOfExtractionFiles > 1)
\t\t\t{
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfileName2 = paste(outputName,"_LR-GLM-CA_results.txt",sep="")
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tfileName2 = paste(outputName,"_LR-GLM_results.txt",sep="")
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tfileName2 = paste(outputName, "_linear_regression_results.txt", sep="")
\t\t\t\t\t}
\t\t\t\tuniLRcoefficientsNames1 = c(); uniLRRsquaresNames1 = c(); uniLRdeltaRsquaresNames1 = c(); uniLRpValuesNames = c()
\t\t\t\tuniLRcoefficientsNames2 = c(); uniLRRsquaresNames2 = c(); uniLRdeltaRsquaresNames2 = c()
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMcoefficientsNames = c(); multiGLMpValuesNames = c()
\t\t\t\t\t\tmat = cbind(realUniLRcoefficients1[,1:dim(realUniLRcoefficients1)[2]],realUniLRRsquares1[,1:dim(realUniLRRsquares1)[2]],
\t\t\t\t\t\t\t\t\trealUniDeltaRsquares1[,2:dim(realUniDeltaRsquares1)[2]],realMultiGLMcoefficients[,1:dim(realMultiGLMcoefficients)[2]])
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsNames = c(); multiGLMCAcommonContributionsNames = c()
\t\t\t\t\t\t\t\tmat = cbind(mat, realMultiGLMCAuniqueContributions[,1:dim(realMultiGLMCAuniqueContributions)[2]],
\t\t\t\t\t\t\t\t\t\t\t\t realMultiGLMCAcommonContributions[,1:dim(realMultiGLMCAcommonContributions)[2]])
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tmat = cbind(realUniLRcoefficients1[,1:dim(realUniLRcoefficients1)[2]], realUniLRRsquares1[,1:dim(realUniLRRsquares1)[2]], 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   realUniDeltaRsquares1[,2:dim(realUniDeltaRsquares1)[2]])
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmat = cbind(mat, realUniLRcoefficients2[,1:dim(realUniLRcoefficients2)[2]], realUniLRRsquares2[,1:dim(realUniLRRsquares2)[2]], 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   realUniDeltaRsquares2[,2:dim(realUniDeltaRsquares2)[2]])
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tuniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("Univariate_LR1_coefficients_", names(envVariables[[h]]), sep=""))\t
\t\t\t\t\t\t\t\tuniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("Univariate_LR1_R2_", names(envVariables[[h]]), sep=""))\t\t\t\t
\t\t\t\t\t\t\t\tuniLRcoefficientsNames2 = cbind(uniLRcoefficientsNames2, paste("Univariate_LR2_coefficients_", names(envVariables[[h]]), sep=""))\t
\t\t\t\t\t\t\t\tuniLRRsquaresNames2 = cbind(uniLRRsquaresNames2, paste("Univariate_LR2_R2_", names(envVariables[[h]]), sep=""))\t
\t\t\t\t\t\t\t\tif (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1, paste("Univariate_LR1_delta_R2_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t\tif (h > 1) uniLRdeltaRsquaresNames2 = cbind(uniLRdeltaRsquaresNames2, paste("Univariate_LR2_delta_R2_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tuniLRcoefficientsNames1 = cbind(uniLRcoefficientsNames1, paste("Univariate_LR_coefficients_", names(envVariables[[h]]), sep=""))\t
\t\t\t\t\t\t\t\tuniLRRsquaresNames1 = cbind(uniLRRsquaresNames1, paste("Univariate_LR_R2_", names(envVariables[[h]]), sep=""))\t\t\t\t
\t\t\t\t\t\t\t\tif (h > 1) uniLRdeltaRsquaresNames1 = cbind(uniLRdeltaRsquaresNames1, paste("Univariate_LR_delta_R2_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMcoefficientsNames = cbind(multiGLMcoefficientsNames, paste("Multivariate_GLM_coefficients_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsNames = cbind(multiGLMCAuniqueContributionsNames,
\t\t\t\t\t\t\t\t\t\t\tpaste0("Multivariate_GLM-CA_unique_contributions_",names(envVariables[[h]])))
\t\t\t\t\t\t\t\t\t\tmultiGLMCAcommonContributionsNames = cbind(multiGLMCAcommonContributionsNames,
\t\t\t\t\t\t\t\t\t\t\tpaste0("Multivariate_GLM-CA_common_contributions_",names(envVariables[[h]])))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tnames = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRpValuesNames, uniLRdeltaRsquaresNames1, multiGLMcoefficientsNames, multiGLMpValuesNames)
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tnames = cbind(names, multiGLMCAuniqueContributionsNames, multiGLMCAcommonContributionsNames)
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tnames = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1, 
\t\t\t\t\t\t\t\t\t\t\t  uniLRcoefficientsNames2, uniLRRsquaresNames2, uniLRdeltaRsquaresNames2)
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tnames = cbind(uniLRcoefficientsNames1, uniLRRsquaresNames1, uniLRdeltaRsquaresNames1)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tcolnames(mat) = names; write.table(mat, file=fileName2, row.names=F, quote=F, sep="\\t")
\t\t\t}
\t}

## 3. randomisation step (permutations of least-cost resistances, torus translation or simulation of environmental rasters)

if (nberOfRandomisations > 0)
\t{\t
\t\tif (impactOnVelocity == TRUE)
\t\t\t{
\t\t\t\tuniLRRsquaresLower1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquaresHigher1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquaresRandomisationPValues1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresLower1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresHigher1 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresRandomisationPValues1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresRandomisationBFs1 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)\t
\t\t\t\tuniLRRsquarePValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquarePValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquarePValuesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquaresLower2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquaresHigher2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRRsquaresRandomisationPValues2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresLower2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresHigher2 = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresRandomisationPValues2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tuniLRdeltaRsquaresRandomisationBFs2 = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)\t
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMcoefficientsLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMcoefficientsHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMcoefficientsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMcoefficientPValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMcoefficientPValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMcoefficientPValuesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAuniqueContributionsLower = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAuniqueContributionsHigher = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAuniqueContributionsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAcommonContributionsLower = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAcommonContributionsHigher = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMCAcommonContributionsRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t}
\t\t\t\t# Matrices for e.g. the "median" p-value tests (number of lines = number of randomisations):\t\t
\t\t\t\tuniLRmedianRsquaresSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\tuniLRmedianRsquarePValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\tuniLRmedianDeltaRsquaresSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\t# uniLRdeltaRsquaresPValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMmedianCoefficientsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsSim = matrix(nrow=nberOfRandomisations, ncol=length(envVariables))
\t\t\t\t\t}
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tfor (i in 1:dim(realMultiGLMcoefficients)[2])
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfor (j in 1:dim(realMultiGLMcoefficients)[1])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (is.na(realMultiGLMcoefficients[j,i]))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsLower[j,i] = NA; multiGLMcoefficientsHigher[j,i] = NA
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsRandomisationPValues[j,i] = NA
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesLower[j,i] = NA; multiGLMcoefficientPValuesHigher[j,i] = NA; 
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesRandomisationPValues[j,i] = NA
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}
\t\tif (impactOnDirection == TRUE)
\t\t\t{
\t\t\t\tmeanEnvValuesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanEnvValuesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanEnvValuesRandomisationPValues = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanEnvValuesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
\t\t\t\trateOfPositiveDifferencesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trateOfPositiveDifferencesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trateOfPositiveDifferencesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\trateOfPositiveDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)
\t\t\t\tmeanDifferencesLower = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanDifferencesHigher = matrix(0, nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanDifferencesRandomisationPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\tmeanDifferencesRandomisationBFs = matrix(nrow=length(envVariables), ncol=nberOfRandomisations)\t
\t\t\t}
\t\tfor (s in 1:nberOfRandomisations)
\t\t\t{
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tuniLRdeltaRsquaresRand = matrix(0, nrow=nberOfExtractionFiles, ncol=nberOfRandomisations)
\t\t\t\t\t\tdistancesSim = list()
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tdistancesSim[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\t# 3.1. randomisation step:
\t\t\t\tif ((externalRandomisations == TRUE)|(externalSimulations == TRUE))
\t\t\t\t\t{
\t\t\t\t\t\tsimRasters = list(); simRasters = hullRasters
\t\t\t\t\t\tif (externalRandomisations == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tcat("Analysis on external randomisation ",s,"\\n",sep="")
\t\t\t\t\t\t\t\textractionFileName = "TreeRandomisation"
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (externalSimulations == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tcat("Analysis on external simulation ",s,"\\n",sep="")
\t\t\t\t\t\t\t\textractionFileName = "TreeSimulations"
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (nchar(localTreesDirectory) == 0)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tdata = read.csv(paste(extractionFileName,"_1.csv",sep=""), header=T, dec=".")
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tdata = read.csv(paste(localTreesDirectory,"/",extractionFileName,"_1.csv",sep=""), header=T, dec=".")\t
\t\t\t\t\t\t\t}
\t\t\t\t\t\tdata = data[with(data, order(startYear,endYear)),]
\t\t\t\t\t\tnode1 = list()
\t\t\t\t\t\tnode2 = list()
\t\t\t\t\t\tstartYear = list()
\t\t\t\t\t\tdispersalTime = list()
\t\t\t\t\t\ttreeIDs = list()
\t\t\t\t\t\tdispersalRate = list()
\t\t\t\t\t\tfromCoor = list()
\t\t\t\t\t\ttoCoor = list()
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (t != 1)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (nchar(localTreesDirectory) == 0)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tfileName = paste(extractionFileName,"_",t,".csv",sep="")
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tfileName = paste(localTreesDirectory,"/",extractionFileName,"_",t,".csv", sep="")
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\tdata = read.csv(fileName, h = T)
\t\t\t\t\t\t\t\t\t\tdata = data[with(data, order(endYear, startYear)),]
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tancestralNodeNAonNullRaster = TRUE
\t\t\t\t\t\t\t\twhile (ancestralNodeNAonNullRaster == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tancestralNodeNAonNullRaster = FALSE
\t\t\t\t\t\t\t\t\t\tancestralBranches = which(!data[,"node1"]%in%data[,"node2"])
\t\t\t\t\t\t\t\t\t\tindicesOfBranchesToRemove = c()
\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(ancestralBranches))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(nullRaster, cbind(data[ancestralBranches[i],"startLon"],data[ancestralBranches[i],"startLat"]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeNAonNullRaster = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\tindicesOfBranchesToRemove = c(indicesOfBranchesToRemove, ancestralBranches[i])
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (length(indicesOfBranchesToRemove) > 0)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tdata = data[-indicesOfBranchesToRemove,]
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (onlyTipBranches == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tindices = which(!data[,"node2"]%in%data[,"node1"])
\t\t\t\t\t\t\t\t\t\tdata = data[indices,]
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tnberOfConnections[t] = dim(data)[1]\t\t
\t\t\t\t\t\t\t\tnode1[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\t\t\t\t\t\tnode1[[t]][] = data[,"node1"]
\t\t\t\t\t\t\t\tnode2[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\t\t\t\t\t\tnode2[[t]][] = data[,"node2"]
\t\t\t\t\t\t\t\tstartYear[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)\t
\t\t\t\t\t\t\t\tstartYear[[t]][] = data[,"startYear"]
\t\t\t\t\t\t\t\tdispersalTime[[t]] = matrix(nrow=nberOfConnections[t], ncol=1)
\t\t\t\t\t\t\t\tdispersalTime[[t]][] = (data[,"endYear"]-data[,"startYear"])
\t\t\t\t\t\t\t\tcolnames(dispersalTime[[t]]) = "dispersalTime"
\t\t\t\t\t\t\t\tfromCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\t\t\t\t\t\t\tfromCoor[[t]][] = cbind(data[,"startLon"], data[,"startLat"])
\t\t\t\t\t\t\t\ttoCoor[[t]] = matrix(nrow=nberOfConnections[t], ncol=2)
\t\t\t\t\t\t\t\ttoCoor[[t]][] = cbind(data[,"endLon"], data[,"endLat"])
\t\t\t\t\t\t\t\ttotalnberOfConnections = totalnberOfConnections + nberOfConnections[t] 
\t\t\t\t\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tdistances[[t]] = matrix(nrow=nberOfConnections[t], ncol=length(envVariables))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (("treeID"%in%colnames(data)) == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttreeIDs[[t]] = data[1,"treeID"]
\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\ttreeIDs[[t]] = "noTreeID"
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor\t
\t\t\t\t\t}
\t\t\t\tif (torusRandomisations == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimRasters = list()
\t\t\t\t\t\tcat("Analysis on torus translation randomisation ",s,"\\n",sep="")\t
\t\t\t\t\t\tfor (h in 2:length(hullRasters))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tsimRasters[[h]] = torusRandomisation(hullRasters[[h]])
\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tplotRaster(merge(simRasters[[h]],envVariables[[h]]), addLegend=T)
\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\ttext1 = paste("torus randomisation ",s," of ",names(hullRasters[[h]])[1],sep="")
\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor\t\t\t\t\t\t
\t\t\t\t\t}
\t\t\t\tif (rastersSimulations == TRUE)
\t\t\t\t\t{\t
\t\t\t\t\t\tsimRasters = list()
\t\t\t\t\t\tcat("Analysis on raster simulation ", s, "\\n", sep="")
\t\t\t\t\t\tfor (h in 2:length(hullRasters))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tsimRasters[[h]] = rasterSimulation(hullRasters[[h]], variogramModels[[h-1]]) 
\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tplotRaster(merge(simRasters[[h]],envVariables[[h]]), addLegend=T)
\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\ttext1 = paste("raster simulation ",s," for ",names(hullRasters[[h]])[1],sep="")
\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t}
\t\t\t\tif (branchRandomisation3 == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimRasters = list()
\t\t\t\t\t\tsimRasters = hullRasters
\t\t\t\t\t\tcat("Analysis of randomised branch positions ",s,"\\n",sep="")
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tcounter1 = 0
\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE; counter1 = counter1+1
\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][,] = NA; toCoorRand[[t]][,] = NA
\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
\t\t\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t\t\t\t# for (i in 1:dim(fromCoor[[t]])[1])
\t\t\t\t\t\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t\t\t\t\t\t# segments(fromCoor[[t]][i,1], fromCoor[[t]][i,2], toCoor[[t]][i,1], toCoor[[t]][i,2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tancestralIndex = list(); ancestralNodes = list(); counter = 0
\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(node1[[t]]))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeBoolean = TRUE
\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(node2[[t]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (node1[[t]][i,1] == node2[[t]][j,1])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeBoolean = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeBoolean == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1
\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralIndex[[counter]] = i
\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodes[[counter]] = node1[[t]][i,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(ancestralIndex))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][ancestralIndex[[i]],1] = fromCoor[[t]][ancestralIndex[[i]],1]
\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][ancestralIndex[[i]],2] = fromCoor[[t]][ancestralIndex[[i]],2]
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tancestralNodes = unique(ancestralNodes)
\t\t\t\t\t\t\t\t\t\tstartingNodes = list(); startingNodes = ancestralNodes # startingNodes[[1]] = ancestralNode[1]
\t\t\t\t\t\t\t\t\t\twhile (length(startingNodes) > 0)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tnewStartingNodes = list(); c = 0
\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(startingNodes))
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tnodes2 = node2[[t]][which(node1[[t]][,1]==startingNodes[[i]]),1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (length(nodes2) > 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(nodes2))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewStartingNodes[[c]] = nodes2[j]\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tk = which(node2[[t]][,1]==nodes2[j])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt01 = c(fromCoor[[t]][k,1], fromCoor[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt02 = c(toCoor[[t]][k,1], toCoor[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = pt02[1]-pt01[1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = pt02[2]-pt01[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoorRand[[t]][k,1], fromCoorRand[[t]][k,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(NA,NA)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[1] = pt1[1]+xTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[2] = pt1[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter2 = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter2 = counter2+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter2 > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print("counter2 > 100")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeID = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(ancestralNodes))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodes[[h]] == node1[[t]][which(node2[[t]][,1]==nodes2[j]),1])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tancestralNodeID = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeID == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter1 <= 10)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (ancestralNodeID == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(c(ancestralNodeID,startingNodes[[i]],nodes2[j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2[1], pt2[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][k,1] = pt1[1]; fromCoorRand[[t]][k,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][k,1] = pt2[1]; toCoorRand[[t]][k,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttoModify = which(node1[[t]][,1]==nodes2[j])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (k in 1:length(toModify))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][toModify[k],1] = pt2[1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][toModify[k],2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\tstartingNodes = newStartingNodes\t
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (branchRandomisation2 == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimRasters = list()
\t\t\t\t\t\tsimRasters = hullRasters
\t\t\t\t\t\tcat("Analysis of randomised branch positions ", s, "\\n", sep="")
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=F)
\t\t\t\t\t\t\t\t\t\tlines(points[hull,], lwd=0.5, col="black")
\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tfor (i in 1:nberOfConnections[t])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (rotatingEndNodes == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\tcounter = 0
\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE; twoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1[1], pt1[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (rotatingEndNodes == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt1[1]; toCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt2[1]; fromCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}\t\t\t
\t\t\t\tif (branchRandomisation1 == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimRasters = list()
\t\t\t\t\t\tsimRasters = hullRasters
\t\t\t\t\t\tcat("Analysis of randomised branch positions ", s, "\\n", sep="")
\t\t\t\t\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (t == 1) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\tif (t >= 2) plotRaster(envVariables[[1]], addLegend=T, new=T)
\t\t\t\t\t\t\t\t\t\tplot(sps, lwd=0.5, border="black", add=T)
\t\t\t\t\t\t\t\t\t\ttext1 = paste("randomisation of branch positions, sampled tree ", t, sep="")
\t\t\t\t\t\t\t\t\t\tmtext(text1, col="black", cex=0.7, line=0)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tfor (i in 1:nberOfConnections[t])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\twhile (twoPointsOnTheGrid == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tpt1NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\twhile (pt1NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated = pt1
\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = runif(1)*(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = runif(1)*(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[1] = pt1[1]+xTranslation; pt1_translated[2] = pt1[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt1_translated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[1] = pt1_translated[1]-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txTranslation = xTranslation-(extent(hullRasters[[1]])@xmax-extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt1_translated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1_translated[2] = pt1_translated[2]-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyTranslation = yTranslation-(extent(hullRasters[[1]])@ymax-extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt1_translated[1],pt1_translated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tinsideAtLeastOneHullPolygon = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(sps))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[h]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt1_translated[1], pt1_translated[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tinsideAtLeastOneHullPolygon = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (insideAtLeastOneHullPolygon == FALSE) NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = pt1_translated
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2[1] = pt2[1]+xTranslation; pt2[2] = pt2[2]+yTranslation
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tpol_index = NA
\t\t\t\t\t\t\t\t\t\t\t\tfor (j in 1:length(sps))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[j]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt1[1], pt1[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_index = j
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\tcounter = 0
\t\t\t\t\t\t\t\t\t\t\t\twhile (pt2NAarea == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tcounter = counter+1\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(hullRasters[[1]])@xmax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(hullRasters[[1]])@xmin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(hullRasters[[1]])@ymax)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(hullRasters[[1]])@ymin)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(hullRasters))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(hullRasters[[h]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_x = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpol_y = sps[pol_index]@polygons[[1]]@Polygons[[1]]@coords[,2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (point.in.polygon(pt2_rotated[1], pt2_rotated[2], pol_x, pol_y) == 1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (counter > 100)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt1 = c(fromCoor[[t]][i,1],fromCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2 = c(toCoor[[t]][i,1],toCoor[[t]][i,2])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tpt2 = pt2_rotated
\t\t\t\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt1_translated[1], pt1_translated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\tpoints(pt2_rotated[1], pt2_rotated[2], pch=16, col="black", cex=0.25)
\t\t\t\t\t\t\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col="black", lwd=0.3)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][i,1] = pt1[1]; fromCoorRand[[t]][i,2] = pt1[2]
\t\t\t\t\t\t\t\t\t\ttoCoorRand[[t]][i,1] = pt2[1]; toCoorRand[[t]][i,2] = pt2[2]
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}\t\t
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tif (distPermutations == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tcat("Analysis on environmental distances permutation ", s, "\\n", sep="")
\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tdistancesSim[[t]][,h] = distances[[t]][sample(dim(distances[[t]])[1]),h]
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tdistancesSim[[t]][,1] = distances[[t]][,1]
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (fourCells == TRUE) directions = 4
\t\t\t\t\t\t\t\t\t\tif (fourCells == FALSE) directions = 8
\t\t\t\t\t\t\t\t\t\tif (resistances[h] == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="c", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], mean, directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="c", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariable = transition(simRasters[[h]], function(x) 1/mean(x), directions)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimTrEnvVariableCorr = geoCorrection(simTrEnvVariable, type="r", multpl=F, scl=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tbuffer = list()
\t\t\t\t\t\t\t\t\t\tbuffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {
\t\t\t\t\t\t\t\t\t\t# for (t in 1:nberOfExtractionFiles) {
\t\t\t\t\t\t\t\t\t\t\t\tmat = matrix(nrow=nberOfConnections[t], ncol=1)
\t\t\t\t\t\t\t\t\t\t\t\tif (straightLineDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlinesList = list()
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoorRand[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpoints = rbind(fromCoorRand[[t]][i,], toCoorRand[[t]][i,])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tlinesList[[i]] = Lines(list(Line(points)),i)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tlines = SpatialLines(linesList)
\t\t\t\t\t\t\t\t\t\t\t\t\t\textractions = raster::extract(simRasters[[h]], lines)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoorRand[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[i] = sum(extractions[[i]], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\tif (leastCostDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[] = diag(costDistance(simTrEnvVariableCorr, fromCoorRand[[t]], toCoorRand[[t]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (randomWalkDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tbranchesNotNA = which(!((is.na(raster::extract(simRasters[[h]],fromCoorRand[[t]][])))|
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t   (is.na(raster::extract(simRasters[[h]],toCoorRand[[t]][])))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimRasterName = paste("CS_rasters/",names(simRasters[[h]]),"_",outputName,"_cs",extensions[h],sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (juliaCSImplementation == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[branchesNotNA,] = circuitScape1(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], fourCells, 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][branchesNotNA,], toCoorRand[[t]][branchesNotNA,], OS, outputName,
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tt, nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[branchesNotNA,] = circuitScape2(simRasters[[h]], simRasterName, resistances[[h]], avgResistances[[h]], fourCells, 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfromCoorRand[[t]][branchesNotNA,], toCoorRand[[t]][branchesNotNA,], OS, outputName,
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tt, nberOfCores_CS)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (commuteDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(fromCoorRand[[t]][,1]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tspatialPoints = SpatialPoints(cbind(c(fromCoorRand[[t]][i,1],toCoorRand[[t]][i,1]),
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tc(fromCoorRand[[t]][i,2],toCoorRand[[t]][i,2])))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[i] = commuteDistance(simTrEnvVariableCorr, spatialPoints)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (rSPDistance == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmat[] = diag(rSPDistance(simTrEnvVariableCorr, fromCoorRand[[t]], toCoorRand[[t]],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t theta=thetaValue, totalNet="total", method=1))
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\t\tcolnames(mat) = names(envVariables[[h]])
\t\t\t\t\t\t\t\t\t\t\t\t# buffer[[t]] = mat
\t\t\t\t\t\t\t\t\t\t\t\tmat
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\tfor (t in 1:length(buffer))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tbuffer[[t]][!is.finite(buffer[[t]][])] = NA
\t\t\t\t\t\t\t\t\t\t\t\tbuffer[[t]][buffer[[t]][]==-1] = NA
\t\t\t\t\t\t\t\t\t\t\t\tdistancesSim[[t]][,h] = buffer[[t]][]
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (impactOnDirection == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimMeanExtractions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimRateOfPositiveDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimMeanDifferences = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tenvValues = 0; ancestralNodes = unique(node1[[t]][which(!node1[[t]]%in%node2[[t]])])
\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(ancestralNodes))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tancestralBranch = which(node1[[t]]==ancestralNodes[i])[1]
\t\t\t\t\t\t\t\t\t\t\t\tenvValues = envValues + raster::extract(envVariables[[h]],
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcbind(fromCoorRand[[t]][ancestralBranch,1],fromCoorRand[[t]][ancestralBranch,2]))
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tsimMeanExtractions[t,h] = envValues + mean(raster::extract(envVariables[[h]], cbind(toCoorRand[[t]][,1],toCoorRand[[t]][,2])), na.rm=T)
\t\t\t\t\t\t\t\t\t\tif ((!is.na(meanEnvValues[t,h])) & (!is.na(simMeanExtractions[t,h])))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (meanEnvValues[t,h] < simMeanExtractions[t,h]) meanEnvValuesHigher[t,h] = meanEnvValuesHigher[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\tif (meanEnvValues[t,h] > simMeanExtractions[t,h]) meanEnvValuesLower[t,h] = meanEnvValuesLower[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tmeanEnvValuesHigher[t,h] = NA; meanEnvValuesLower[t,h] = NA
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tdiffs = raster::extract(envVariables[[h]], fromCoorRand[[t]])-raster::extract(envVariables[[h]], toCoorRand[[t]])
\t\t\t\t\t\t\t\t\t\tsimRateOfPositiveDifferences[t,h] = sum(diffs[!is.na(diffs)] > 0)/length(diffs[!is.na(diffs)])
\t\t\t\t\t\t\t\t\t\tsimMeanDifferences[t,h] = mean(diffs, na.rm=T)
\t\t\t\t\t\t\t\t\t\tif ((!is.na(rateOfPositiveDifferences[t,h])) & (!is.na(sum(diffs > 0))))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (rateOfPositiveDifferences[t,h] < simRateOfPositiveDifferences[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\trateOfPositiveDifferencesHigher[t,h] = rateOfPositiveDifferencesHigher[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (rateOfPositiveDifferences[t,h] > simRateOfPositiveDifferences[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\trateOfPositiveDifferencesLower[t,h] = rateOfPositiveDifferencesLower[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\trateOfPositiveDifferencesLower[t,h] = NA; rateOfPositiveDifferencesHigher[t,h] = NA
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif ((!is.na(meanDifferences[t,h])) & (!is.na(mean(diffs))))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (meanDifferences[t,h] < simMeanDifferences[t,h]) meanDifferencesHigher[t,h] = meanDifferencesHigher[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\tif (meanDifferences[t,h] > simMeanDifferences[t,h]) meanDifferencesLower[t,h] = meanDifferencesLower[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmeanDifferencesLower[t,h] = NA; meanDifferencesHigher[t,h] = NA
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\t
\t\t\t\t# 3.2. Linear regressions analyses on randomised datasets:\t
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tsimUniLRRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimUniLRRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimUniLRRsquarePValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimUniLRDeltaRsquares1 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tsimUniLRDeltaRsquares2 = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tsimMultiGLMcoefficients = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\t\t\tsimMultiGLMcoefficientPValues = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\t\t\tsimMultiGLMCAuniqueContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\t\t\tsimMultiGLMCAcommonContributions = matrix(nrow=nberOfExtractionFiles, ncol=length(envVariables))
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (distPermutations == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tdistVariables = paste("dispersalTime[[",t,"]]"," ~ distancesSim[[",t,"]][,1]",sep="")
\t\t\t\t\t\t\t\t\t\tLM = lm(as.formula(distVariables))
\t\t\t\t\t\t\t\t\t\tsimUniLRRsquares1[t,1] = summary(LM)$r.squared
\t\t\t\t\t\t\t\t\t\tf = summary(LM)$fstatistic
\t\t\t\t\t\t\t\t\t\tif (is.numeric(f))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tp = pf(f[1],f[2],f[3],lower.tail=F)
\t\t\t\t\t\t\t\t\t\t\t\tattributes(p) = NULL; simUniLRRsquarePValues[t,h] = p
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tsimUniLRRsquarePValues[t,h] = NA 
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tdistVariables = paste("dispersalTime[[",t,"]]"," ~ distancesSim[[",t,"]][,",h,"]",sep="")
\t\t\t\t\t\t\t\t\t\tLM = lm(as.formula(distVariables))
\t\t\t\t\t\t\t\t\t\tsimUniLRRsquares1[t,h] = summary(LM)$r.squared
\t\t\t\t\t\t\t\t\t\tf = summary(LM)$fstatistic; p = pf(f[1],f[2],f[3],lower.tail=F)
\t\t\t\t\t\t\t\t\t\tattributes(p) = NULL; simUniLRRsquarePValues[t,h] = p
\t\t\t\t\t\t\t\t\t\tsimUniLRDeltaRsquares1[t,h] = simUniLRRsquares1[t,h]-simUniLRRsquares1[t,1]
\t\t\t\t\t\t\t\t\t\tnonNaNdispersalTimes = dispersalTime[[t]][!is.na(distancesSim[[t]][,1])]
\t\t\t\t\t\t\t\t\t\tnonNaNdistancesSim1 = distancesSim[[t]][,1][!is.na(distancesSim[[t]][,1])]
\t\t\t\t\t\t\t\t\t\tnonNaNdistancesSim2 = distancesSim[[t]][,h][!is.na(distancesSim[[t]][,1])]
\t\t\t\t\t\t\t\t\t\tdistVariables = paste("nonNaNdistancesSim2 ~ nonNaNdistancesSim1",sep="")
\t\t\t\t\t\t\t\t\t\tLM1 = lm(as.formula(distVariables)); residuals_LM1 = LM1$residuals; indices = which(!is.na(residuals_LM1))
\t\t\t\t\t\t\t\t\t\tdistVariables = paste("nonNaNdispersalTimes[indices] ~ nonNaNdistancesSim2[indices] + residuals_LM1",sep="")
\t\t\t\t\t\t\t\t\t\tLM2 = lm(as.formula(distVariables)); simUniLRRsquares2[t,h] = summary(LM2)$r.squared
\t\t\t\t\t\t\t\t\t\tsimUniLRDeltaRsquares2[t,h] = simUniLRRsquares2[t,h]-simUniLRRsquares1[t,1]
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (all == FALSE)
\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\tfor (i in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistVariables = paste("dispersalTime ~ ", names(envVariables[[1]]), sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,1] = dispersalTime[[t]]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,2] = distances[[t]][,1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c("dispersalTime")
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c(matMultiGLMNames, names(envVariables[[1]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (h == i)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,h+1] = distancesSim[[t]][,h]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,h+1] = distances[[t]][,h]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\t\t\t\t\t\t\t\t\tm = min(unlist(lapply(as.matrix(matMultiGLM),min)), na.rm=T)\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, preLogTransformation, m)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, logTransformation)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, zTransformation)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\t\t\t\t\t\t\t\t\tform = as.formula(distVariables)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLM = stats::glm(form, data=matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tCA = calc.yhat(multiGLM, prec=5)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnames(multiGLM$coefficients)[h] = names(envVariables[[h]])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMcoefficients[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMcoefficientPValues[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMCAuniqueContributions[t,i] = CA$PredictorMetrics[h,"Unique"]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMCAcommonContributions[t,i] = CA$PredictorMetrics[h,"Common"]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tdistVariables = paste("dispersalTime ~ ", names(envVariables[[1]]), sep="")
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = matrix(nrow=length(dispersalTime[[t]]), ncol=(length(envVariables)+1))
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,1] = dispersalTime[[t]]
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,2] = distances[[t]][,1]
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c("dispersalTime")
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c(matMultiGLMNames, names(envVariables[[1]]))
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLMNames = c(matMultiGLMNames, names(envVariables[[h]]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\tdistVariables = paste(distVariables, " + ", names(envVariables[[h]]), sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\t# distVariables = paste(distVariables, " + ", "distancesSim[[", t, "]][,", h, "]", sep="")
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM[,h+1] = distancesSim[[t]][,h]\t
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = lapply(matMultiGLM, zTransformation)
\t\t\t\t\t\t\t\t\t\t\t\tmatMultiGLM = as.data.frame(matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\tnames(matMultiGLM) = matMultiGLMNames
\t\t\t\t\t\t\t\t\t\t\t\tform = as.formula(distVariables)
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLM = stats::glm(form, data=matMultiGLM)
\t\t\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tCA = calc.yhat(multiGLM, prec=5)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tnames(multiGLM$coefficients)[h] = names(envVariables[[h]])
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMcoefficients[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Estimate"]
\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMcoefficientPValues[t,i] = summary(multiGLM)[]$coefficients[names(envVariables[[h]]),"Pr(>|t|)"]
\t\t\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMCAuniqueContributions[t,i] = CA$PredictorMetrics[h,"Unique"]
\t\t\t\t\t\t\t\t\t\t\t\t\t\tsimMultiGLMCAcommonContributions[t,i] = CA$PredictorMetrics[h,"Common"]
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (distPermutations == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tH = 1
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tH = 2
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (impactOnDirection == TRUE) H = 2

\t\t\t\t# 3.3.1. Count of lower & higher values compared to real/observed values for this randomisation:\t\t\t
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (simUniLRRsquares1[t,h] < realUniLRRsquares1[t,h])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquaresLower1[t,h] = uniLRRsquaresLower1[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquaresHigher1[t,h] = uniLRRsquaresHigher1[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (simUniLRRsquarePValues[t,h] < realUniLRRsquarePValues[t,h])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquarePValuesLower[t,h] = uniLRRsquarePValuesLower[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquarePValuesHigher[t,h] = uniLRRsquarePValuesHigher[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (simUniLRDeltaRsquares1[t,h] < realUniDeltaRsquares1[t,h])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresLower1[t,h] = uniLRdeltaRsquaresLower1[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresHigher1[t,h] = uniLRdeltaRsquaresHigher1[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (simUniLRRsquares2[t,h] < realUniLRRsquares2[t,h])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquaresLower2[t,h] = uniLRRsquaresLower2[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRRsquaresHigher2[t,h] = uniLRRsquaresHigher2[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (simUniLRDeltaRsquares2[t,h] < realUniDeltaRsquares2[t,h])
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresLower2[t,h] = uniLRdeltaRsquaresLower2[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresHigher2[t,h] = uniLRdeltaRsquaresHigher2[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (!is.na(simMultiGLMcoefficients[t,h]) & !is.na(realMultiGLMcoefficients[t,h]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMultiGLMcoefficients[t,h] < realMultiGLMcoefficients[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsLower[t,h] = multiGLMcoefficientsLower[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsHigher[t,h] = multiGLMcoefficientsHigher[t,h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\t\t\tif (!is.na(simMultiGLMcoefficientPValues[t,h]) & !is.na(realMultiGLMcoefficientPValues[t,h]))
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMultiGLMcoefficientPValues[t,h] < realMultiGLMcoefficientPValues[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesLower[t,h] = multiGLMcoefficientPValuesLower[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesHigher[t,h] = multiGLMcoefficientPValuesHigher[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMultiGLMCAcommonContributions[t,h] < realMultiGLMCAcommonContributions[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMCAcommonContributionsLower[t,h] = multiGLMCAcommonContributionsLower[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMCAcommonContributionsHigher[t,h] = multiGLMCAcommonContributionsHigher[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMultiGLMCAuniqueContributions[t,h] < realMultiGLMCAuniqueContributions[t,h])
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsLower[t,h] = multiGLMCAuniqueContributionsLower[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsHigher[t,h] = multiGLMCAuniqueContributionsHigher[t,h] + 1;
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\t\t
\t\t\t\t# 3.3.2. Computations of the median and BF values for this randomisation:\t
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tuniLRmedianRsquaresSim[s,h] = median(simUniLRRsquares1[,h])
\t\t\t\t\t\t\t\tuniLRmedianRsquarePValuesSim[s,h] = median(simUniLRRsquarePValues[,h])
\t\t\t\t\t\t\t\tuniLRmedianDeltaRsquaresSim[s,h] = median(simUniLRRsquares1[,h]-simUniLRRsquares1[,1])
\t\t\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientsSim[s,h] = median(simMultiGLMcoefficients[,h])
\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesSim[s,h] = median(simMultiGLMcoefficientPValues[,h])
\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsSim[s,h] = median(simMultiGLMCAuniqueContributions[,h])
\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsSim[s,h] = median(simMultiGLMCAcommonContributions[,h])
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tc = 0
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (simUniLRDeltaRsquares1[t,h] < realUniDeltaRsquares1[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tf = c/nberOfExtractionFiles; bf = f/(1-f)
\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationBFs1[h,s] = round(bf, 4)
\t\t\t\t\t\t\t\tc = 0
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (simUniLRDeltaRsquares2[t,h] < realUniDeltaRsquares2[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tf = c/nberOfExtractionFiles; bf = f/(1-f)
\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationBFs2[h,s] = round(bf, 4)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (showingPlots == TRUE) dev.off()
\t\t\t\t\t}
\t\t\t\tif (impactOnDirection == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tc = 0; missingValues = 0
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif ((!is.na(simMeanExtractions[t,h]))&(!is.na(meanEnvValues[t,h])))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (resistances[h] == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMeanExtractions[t,h] > meanEnvValues[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMeanExtractions[t,h] < meanEnvValues[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmissingValues = missingValues + 1
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tf = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
\t\t\t\t\t\t\t\tmeanEnvValuesRandomisationBFs[h,s] = round(bf, 4)
\t\t\t\t\t\t\t\tc = 0; missingValues = 0
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif ((!is.na(simMeanDifferences[t,h]))&(!is.na(meanDifferences[t,h])))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (resistances[h] == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMeanDifferences[t,h] < meanDifferences[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simMeanDifferences[t,h] > meanDifferences[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmissingValues = missingValues + 1
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tf = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
\t\t\t\t\t\t\t\tmeanDifferencesRandomisationBFs[h,s] = round(bf, 4)
\t\t\t\t\t\t\t\tc = 0; missingValues = 0
\t\t\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif ((!is.na(simRateOfPositiveDifferences[t,h])) & (!is.na(rateOfPositiveDifferences[t,h])))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tif (resistances[h] == TRUE)
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simRateOfPositiveDifferences[t,h] < rateOfPositiveDifferences[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (simRateOfPositiveDifferences[t,h] > rateOfPositiveDifferences[t,h]) c = c+1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tmissingValues = missingValues + 1
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\tf = c/(nberOfExtractionFiles-missingValues); bf = f/(1-f)
\t\t\t\t\t\t\t\trateOfPositiveDifferencesRandomisationBFs[h,s] = round(bf, 4)
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t} # End of randomisations
\t\t\t\t
\t\t# 3.4.1. Computations of p-values for the histogram constructions:\t
\t\tif (impactOnVelocity == TRUE)
\t\t\t{\t
\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tuniLRRsquaresRandomisationPValues1[t,h] = (nberOfRandomisations-uniLRRsquaresLower1[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\tuniLRRsquarePValuesRandomisationPValues[t,h] = (uniLRRsquarePValuesLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationPValues1[t,h] = (nberOfRandomisations-uniLRdeltaRsquaresLower1[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\tuniLRRsquaresRandomisationPValues2[t,h] = (nberOfRandomisations-uniLRRsquaresLower2[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationPValues2[t,h] = (nberOfRandomisations-uniLRdeltaRsquaresLower2[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMcoefficientsLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesRandomisationPValues[t,h] = (multiGLMcoefficientPValuesLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMCAuniqueContributionsLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\tmultiGLMCAcommonContributionsRandomisationPValues[t,h] = (nberOfRandomisations-multiGLMCAcommonContributionsLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}\t
\t\tif (impactOnDirection == TRUE)
\t\t\t{\t
\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tfor (t in 1:nberOfExtractionFiles)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (resistances[h] == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmeanEnvValuesRandomisationPValues[t,h] = (meanEnvValuesLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\trateOfPositiveDifferencesRandomisationPValues[t,h] = (rateOfPositiveDifferencesHigher[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\tmeanDifferencesRandomisationPValues[t,h] = (meanDifferencesHigher[t,h])/nberOfRandomisations\t
\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\tmeanEnvValuesRandomisationPValues[t,h] = (meanEnvValuesHigher[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\trateOfPositiveDifferencesRandomisationPValues[t,h] = (rateOfPositiveDifferencesLower[t,h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t\tmeanDifferencesRandomisationPValues[t,h] = (meanDifferencesLower[t,h])/nberOfRandomisations\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}
\t\t
\t\t# 3.4.2. Computations of the "median" p-values (one-tailed test, except for multiGLM and CA coefficients: two-tailed tests):\t\t
\t\tif (impactOnVelocity == TRUE)
\t\t\t{
\t\t\t\tuniLRmedianRsquaresLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\tuniLRmedianRsquaresPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tuniLRmedianRsquarePValuesLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\tuniLRmedianRsquarePValuesPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tuniLRmedianDeltaRsquaresLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\tuniLRmedianDeltaRsquaresPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMmedianCoefficientsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))\t
\t\t\t\t\t\tmultiGLMmedianCoefficientsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCoefficientsPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesLowerValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))\t
\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsLowerValues = matrix(0, nrow=1, ncol=length(envVariables))\t
\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsHigherValues = matrix(0, nrow=1, ncol=length(envVariables))
\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsPValues = matrix(nrow=1, ncol=length(envVariables))
\t\t\t\t\t}
\t\t\t\tfor (h in H:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tfor (s in 1:nberOfRandomisations)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (uniLRmedianRsquaresSim[s,h] < uniLRmedianRsquares[h])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tuniLRmedianRsquaresLowerValues[h] = uniLRmedianRsquaresLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (uniLRmedianRsquarePValuesSim[s,h] < uniLRmedianRsquarePValues[h])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tuniLRmedianRsquarePValuesLowerValues[h] = uniLRmedianRsquarePValuesLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (uniLRmedianDeltaRsquaresSim[s,h] < uniLRmedianDeltaRsquares[h])
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tuniLRmedianDeltaRsquaresLowerValues[h] = uniLRmedianDeltaRsquaresLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (!is.na(multiGLMmedianCoefficients[h]))
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tif (multiGLMmedianCoefficientsSim[s,h] < multiGLMmedianCoefficients[h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientsLowerValues[h] = multiGLMmedianCoefficientsLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientsHigherValues[h] = multiGLMmedianCoefficientsHigherValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (multiGLMmedianCoefficientPValuesSim[s,h] < multiGLMmedianCoefficientPValues[h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesLowerValues[h] = multiGLMmedianCoefficientPValuesLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tif (multiGLMmedianCAuniqueContributionsSim[s,h] < multiGLMmedianCAuniqueContributions[h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsLowerValues[h] = multiGLMmedianCAuniqueContributionsLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsHigherValues[h] = multiGLMmedianCAuniqueContributionsHigherValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tif (multiGLMmedianCAcommonContributionsSim[s,h] < multiGLMmedianCAcommonContributions[h])
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsLowerValues[h] = multiGLMmedianCAcommonContributionsLowerValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsHigherValues[h] = multiGLMmedianCAcommonContributionsHigherValues[h] + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}\t\t\t\t
\t\t\t\t\t\t\t}
\t\t\t\t\t\tuniLRmedianRsquaresPValues[h] = (nberOfRandomisations-uniLRmedianRsquaresLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\tuniLRmedianRsquarePValuesPValues[h] = (uniLRmedianRsquarePValuesLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\tuniLRmedianDeltaRsquaresPValues[h] = (nberOfRandomisations-uniLRmedianDeltaRsquaresPValues[h])/nberOfRandomisations
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientsPValues[h] = (nberOfRandomisations-multiGLMmedianCoefficientsLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\t\t\tmultiGLMmedianCoefficientPValuesPValues[h] = (multiGLMmedianCoefficientPValuesLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMmedianCAuniqueContributionsPValues[h] = (nberOfRandomisations-multiGLMmedianCAuniqueContributionsLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\t\t\tmultiGLMmedianCAcommonContributionsPValues[h] = (nberOfRandomisations-multiGLMmedianCAcommonContributionsLowerValues[h])/nberOfRandomisations
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}
\t\t
\t\t# 3.5. Creation of the pdf file containing histograms:
\t\tif (impactOnVelocity == TRUE)
\t\t\t{
\t\t\t\tfileName = paste(outputName, "_randomisation_results.pdf", sep="")\t
\t\t\t\tif ((plottingHistograms == TRUE)&(nberOfExtractionFiles > 49)&(nberOfRandomisations > 1))
\t\t\t\t\t{
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tpdf(fileName, width=(4*length(envVariables)), height=(5*4))
\t\t\t\t\t\t\t\t\t\tpar(mfrow=c(7,length(envVariables)))
\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\tpdf(fileName, width=(4*length(envVariables)), height=(3*4))
\t\t\t\t\t\t\t\t\t\tpar(mfrow=c(5,length(envVariables)))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tpdf(fileName, width=(4*length(envVariables)), height=(2.25*4))
\t\t\t\t\t\t\t\tpar(mfrow=c(3,length(envVariables)))
\t\t\t\t\t\t\t}
\t\t\t\t\t\tbreakList = (0:50)/50
\t\t\t\t\t\txMin = min(realUniLRRsquares1[,1], na.rm=T); xMax = max(realUniLRRsquares1[,1], na.rm=T)\t
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tif (xMin > min(realUniLRRsquares1[,h], na.rm=T)) xMin = min(realUniLRRsquares1[,h], na.rm=T)
\t\t\t\t\t\t\t\tif (xMax < max(realUniLRRsquares1[,h], na.rm=T)) xMax = max(realUniLRRsquares1[,h], na.rm=T)

\t\t\t\t\t\t\t}
\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmin = min(realUniLRRsquares1[,h]); max = max(realUniLRRsquares1[,h])
\t\t\t\t\t\t\t\thist(realUniLRRsquares1[,h], freq=T, xlim=c(0,xMax), breaks=seq(0,xMax,by=(xMax-0)/50), xlab="Univariate LR R2's", 
\t\t\t\t\t\t\t\t\t main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (distPermutations == FALSE)
\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ta = hist(uniLRRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ta = hist(uniLRRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\thist(uniLRRsquaresRandomisationPValues1[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
\t\t\t\t\t\t\t\t\t xlab="Univariate LR R2 p-values (randomisations)", main="geographical distance", cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)\t
\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\thist(uniLRRsquaresRandomisationPValues1[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
\t\t\t\t\t\t\t\t\t xlab="Univariate LR R2 p-values (randomisations)", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\ta = hist(uniLRdeltaRsquaresRandomisationPValues1[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tplot.new()
\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\thist(uniLRdeltaRsquaresRandomisationPValues1[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t xlab="Univariate LR delta R2 (Q) p-values (randomisations)", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\txMin = min(realMultiGLMcoefficients[,1], na.rm=T); xMax = max(realMultiGLMcoefficients[,1], na.rm=T)
\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (xMin > min(realMultiGLMcoefficients[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMin = min(realMultiGLMcoefficients[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (xMax < max(realMultiGLMcoefficients[,h], na.rm=T))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\txMax = max(realMultiGLMcoefficients[,h], na.rm=T)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tif (is.na(mean(realMultiGLMcoefficients[,h])))
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\thist(realMultiGLMcoefficients[,h], freq=T, xlim=c(xMin,xMax), breaks=seq(xMin,xMax,by=(xMax-xMin)/50),
\t\t\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM coefficients", main=names(envVariables[[h]]), cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (distPermutations == FALSE)
\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMcoefficientsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tyMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t} 
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMcoefficientsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\tyMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t} 
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\thist(multiGLMcoefficientsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM coefficient p-values (randomisations)", main="geographical distance",
\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\thist(multiGLMcoefficientsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax),
\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM coefficient p-values (randomisations)", main=names(envVariables[[h]]),
\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\tif (distPermutations == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMCAuniqueContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t} 
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMCAuniqueContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\thist(multiGLMCAuniqueContributionsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main="geographical distance", 
\t\t\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\thist(multiGLMCAuniqueContributionsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main=names(envVariables[[h]]), 
\t\t\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\tif (distPermutations == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMCAcommonContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\tplot.new()
\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\tyMax = 0\t
\t\t\t\t\t\t\t\t\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\ta = hist(multiGLMCAcommonContributionsRandomisationPValues[,h], breaks=breakList, plot=F)
\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (yMax < max(a$counts)) yMax = max(a$counts)
\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\thist(multiGLMCAcommonContributionsRandomisationPValues[,1], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM CA unique contributions p-values (randomisations)", main="geographical distance", 
\t\t\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\thist(multiGLMCAcommonContributionsRandomisationPValues[,h], freq=T, breaks=breakList, xlim=c(0,1), ylim=c(0,yMax), 
\t\t\t\t\t\t\t\t\t\t\t\t\t xlab="Multivariate GLM CA common contributions p-values (randomisations)", main=names(envVariables[[h]]), 
\t\t\t\t\t\t\t\t\t\t\t\t\t cex.main=1.5, cex.axis=1.2, cex=1.2, cex.lab=1.1)
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tdev.off()
\t\t\t\t\t}
\t\t\t}\t\t\t\t
\t\t\t
\t\t# 3.6. Creation of the text files containing randomisation results:
\t\tif ((nberOfRandomisations > 1) & (impactOnVelocity == TRUE))
\t\t\t{
\t\t\t\tfileName = paste(outputName, "_randomisation_results.txt", sep="")
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmat = cbind(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]],
\t\t\t\t\t\t\t\t\tmultiGLMcoefficientsRandomisationPValues[,2:dim(multiGLMcoefficientsRandomisationPValues)[2]])
\t\t\t\t\t\tif (nberOfExtractionFiles == 1)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\t\t  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\t\t  multiGLMcoefficientsRandomisationPValues[,2:dim(multiGLMcoefficientsRandomisationPValues)[2]]))
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmat = cbind(mat, multiGLMCAuniqueContributionsRandomisationPValues[,2:dim(multiGLMCAuniqueContributionsRandomisationPValues)[2]],
\t\t\t\t\t\t\t\t\t\t\t\t multiGLMCAcommonContributionsRandomisationPValues[,2:dim(multiGLMCAcommonContributionsRandomisationPValues)[2]])
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tmat = cbind(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]])
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmat = cbind(mat, uniLRRsquaresRandomisationPValues2[,2:dim(uniLRRsquaresRandomisationPValues2)[2]], 
\t\t\t\t\t\t\t\t\t\t\tuniLRdeltaRsquaresRandomisationPValues2[,2:dim(uniLRdeltaRsquaresRandomisationPValues2)[2]])
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (nberOfExtractionFiles == 1)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\t\t  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]]))
\t\t\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmat = t(c(uniLRRsquaresRandomisationPValues1[,2:dim(uniLRRsquaresRandomisationPValues1)[2]], 
\t\t\t\t\t\t\t\t\t\t  \t\t  uniLRdeltaRsquaresRandomisationPValues1[,2:dim(uniLRdeltaRsquaresRandomisationPValues1)[2]],
\t\t\t\t\t\t\t\t\t\t  \t\t  uniLRRsquaresRandomisationPValues2[,2:dim(uniLRRsquaresRandomisationPValues2)[2]], 
\t\t\t\t\t\t\t\t\t\t\t\t  uniLRdeltaRsquaresRandomisationPValues2[,2:dim(uniLRdeltaRsquaresRandomisationPValues2)[2]]))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tuniLRRsquarePValuesNames1 = c(); uniLRdeltaRsquarePValuesNames1 = c(); uniLRRsquarePValuePValuesNames = c()
\t\t\t\tuniLRRsquarePValuesNames2 = c(); uniLRdeltaRsquarePValuesNames2 = c()
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tmultiGLMcoefficientPValuesNames = c(); multiGLMRsquarePValuePValuesNames = c()
\t\t\t\t\t\tmultiGLMCAuniqueContributionsPValuesNames = c(); multiGLMCAcommonContributionsPValuesNames = c()
\t\t\t\t\t}
\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tuniLRRsquarePValuesNames1 = cbind(uniLRRsquarePValuesNames1, paste("Uni_LR1_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t\tuniLRdeltaRsquarePValuesNames1 = cbind(uniLRdeltaRsquarePValuesNames1, paste("Uni_LR1_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t\tuniLRRsquarePValuesNames2 = cbind(uniLRRsquarePValuesNames2, paste("Uni_LR2_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t\tuniLRdeltaRsquarePValuesNames2 = cbind(uniLRdeltaRsquarePValuesNames2, paste("Uni_LR2_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tuniLRRsquarePValuesNames1 = cbind(uniLRRsquarePValuesNames1, paste("Uni_LR_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t\tuniLRdeltaRsquarePValuesNames1 = cbind(uniLRdeltaRsquarePValuesNames1, paste("Uni_LR_delta_R2_p-values_",names(envVariables[[h]]),sep=""))
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tmultiGLMcoefficientPValuesNames = cbind(multiGLMcoefficientPValuesNames,
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpaste("Multivariate_GLM_coefficient_p-values_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tmultiGLMCAuniqueContributionsPValuesNames = cbind(multiGLMCAuniqueContributionsPValuesNames, 
\t\t\t\t\t\t\t\t\t\tpaste("Multivariate_GLM_CA_unique_contributions_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t\t\t\tmultiGLMCAcommonContributionsPValuesNames = cbind(multiGLMCAcommonContributionsPValuesNames,
\t\t\t\t\t\t\t\t\t\tpaste("Multivariate_GLM_CA_common_contributions_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (GLM == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tnames = cbind(uniLRRsquarePValuesNames, uniLRdeltaRsquarePValuesNames, uniLRRsquarePValuePValuesNames, 
\t\t\t\t\t\t\t\t\t  multiGLMcoefficientPValuesNames, multiGLMRsquarePValuePValuesNames)
\t\t\t\t\t\tif (commonalityAnalysis == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tnames = cbind(names, multiGLMCAuniqueContributionsPValuesNames, multiGLMCAcommonContributionsPValuesNames)
\t\t\t\t\t\t\t}
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tnames = cbind(uniLRRsquarePValuesNames1, uniLRdeltaRsquarePValuesNames1)
\t\t\t\t\t\tif (alternativeQstat == TRUE) names = cbind(names, uniLRRsquarePValuesNames2, uniLRdeltaRsquarePValuesNames2)
\t\t\t\t\t}
\t\t\t\tcolnames(mat) = names; write.table(mat, file=fileName, row.names=F, quote=F, sep="\\t")
\t\t\t}
\t\tif ((nberOfRandomisations > 1) & (impactOnDirection == TRUE))
\t\t\t{
\t\t\t\tfileName = paste(outputName, "_randomisation_results.txt", sep="")
\t\t\t\tmat = cbind(meanEnvValuesRandomisationPValues[,2:dim(meanEnvValuesRandomisationPValues)[2]], 
\t\t\t\t\t\t\trateOfPositiveDifferencesRandomisationPValues[,2:dim(rateOfPositiveDifferencesRandomisationPValues)[2]])
\t\t\t\tif (nberOfExtractionFiles == 1)
\t\t\t\t\t{
\t\t\t\t\t\tmat = t(c(meanEnvValuesRandomisationPValues[,2:dim(meanEnvValuesRandomisationPValues)[2]], 
\t\t\t\t\t\t\t\t  rateOfPositiveDifferencesRandomisationPValues[,2:dim(rateOfPositiveDifferencesRandomisationPValues)[2]]))
\t\t\t\t\t}
\t\t\t\tmeanEnvValuesPValuesNames = c(); rateOfPositiveDifferencesPValuesNames = c()
\t\t\t\tfor (h in 2:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tmeanEnvValuesPValuesNames = cbind(meanEnvValuesPValuesNames, paste("Direction_E_p-values_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t\trateOfPositiveDifferencesPValuesNames = cbind(rateOfPositiveDifferencesPValuesNames, paste("Direction_R_p-values_", names(envVariables[[h]]), sep=""))
\t\t\t\t\t}
\t\t\t\tnames = cbind(meanEnvValuesPValuesNames, rateOfPositiveDifferencesPValuesNames)
\t\t\t\tcolnames(mat) = names; write.table(mat, file=fileName, row.names=F, quote=F, sep="\\t")
\t\t\t}
\t\tif (nberOfExtractionFiles > 49)
\t\t\t{
\t\t\t\tenvVariablesNames = c()
\t\t\t\tfor (h in 1:length(envVariables))
\t\t\t\t\t{
\t\t\t\t\t\tenvVariablesNames = c(envVariablesNames, names(envVariables[[h]]))
\t\t\t\t\t}
\t\t\t\tcolNames1 = c(); colNames2 = c()
\t\t\t\tfor (s in 1:nberOfRandomisations)
\t\t\t\t\t{
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tcolNames1 = c(colNames1, paste("BFs_deltaR2_LR1_randomisation_",s,sep=""))
\t\t\t\t\t\t\t\tcolNames2 = c(colNames2, paste("BFs_deltaR2_LR2_randomisation_",s,sep=""))
\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\tcolNames1 = c(colNames1, paste("BFs_deltaR2_randomisation_",s,sep=""))
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tif (impactOnVelocity == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\trow.names(uniLRdeltaRsquaresRandomisationBFs1) = envVariablesNames
\t\t\t\t\t\tcolnames(uniLRdeltaRsquaresRandomisationBFs1) = colNames1
\t\t\t\t\t\trow.names(uniLRdeltaRsquaresRandomisationBFs2) = envVariablesNames
\t\t\t\t\t\tcolnames(uniLRdeltaRsquaresRandomisationBFs2) = colNames2
\t\t\t\t\t\tfileName = paste(outputName, "_randomisation_Bayes_factors.txt",sep="")
\t\t\t\t\t\tif (alternativeQstat == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tbuffer = cbind(uniLRdeltaRsquaresRandomisationBFs1,uniLRdeltaRsquaresRandomisationBFs2)
\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\tbuffer = uniLRdeltaRsquaresRandomisationBFs1
\t\t\t\t\t\t\t}
\t\t\t\t\t\twrite.table(buffer, fileName, quote=F,sep="\\t")
\t\t\t\t\t}
\t\t\t\tif (impactOnDirection == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tif (pathModel == -1)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\trow.names(rateOfPositiveDifferencesRandomisationBFs) = envVariablesNames
\t\t\t\t\t\t\t\tcolnames(rateOfPositiveDifferencesRandomisationBFs) = colNames
\t\t\t\t\t\t\t\tfileName = paste(outputName,"_positiveDiffs_Bayes_factors.txt",sep="")
\t\t\t\t\t\t\t\twrite.table(rateOfPositiveDifferencesRandomisationBFs, fileName, quote=F, sep="\\t")
\t\t\t\t\t\t\t\trow.names(meanDifferencesRandomisationBFs) = envVariablesNames
\t\t\t\t\t\t\t\tcolnames(meanDifferencesRandomisationBFs) = colNames
\t\t\t\t\t\t\t\tfileName = paste(outputName,"_meanDifferences_Bayes_factors.txt",sep="")
\t\t\t\t\t\t\t\twrite.table(meanDifferencesRandomisationBFs, fileName, quote=F, sep="\\t")
\t\t\t\t\t\t\t}
\t\t\t\t\t\tif (pathModel == 0)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\trow.names(meanEnvValuesRandomisationBFs) = envVariablesNames
\t\t\t\t\t\t\t\tcolnames(meanEnvValuesRandomisationBFs) = colNames
\t\t\t\t\t\t\t\tfileName = paste(outputName,"_direction_E_Bayes_factors.txt",sep="")
\t\t\t\t\t\t\t\twrite.table(meanEnvValuesRandomisationBFs, fileName, quote=F, sep="\\t")
\t\t\t\t\t\t\t\trow.names(rateOfPositiveDifferencesRandomisationBFs) = envVariablesNames
\t\t\t\t\t\t\t\tcolnames(rateOfPositiveDifferencesRandomisationBFs) = colNames
\t\t\t\t\t\t\t\tfileName = paste(outputName,"_direction_R_Bayes_factors.txt",sep="")
\t\t\t\t\t\t\t\twrite.table(rateOfPositiveDifferencesRandomisationBFs, fileName, quote=F, sep="\\t")
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t}
\t}
}
