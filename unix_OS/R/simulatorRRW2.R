simulatorRRW2 <-
function(envVariable=raster(matrix(c(runif(600,5,10),runif(1000,0,5)),nrow=40,ncol=40)),
\t\t\t\t\t\t resistance=TRUE, scalingValue=1, ancestPosition=c(0.4,0.5), birthRate=0.2,
\t\t\t\t\t\t samplingRate=0.2, startingYear=0, samplingWindow=c(10,50), timeSlice=0.1,
\t\t\t\t\t\t timeIntervale=1, showingPlots=FALSE, extractionOfValuesOnMatrix=FALSE) {
\t
\tenvVariable0 = envVariable
\tmat0 = raster::as.matrix(envVariable)
\tenvVariable[!is.na(envVariable[])] = envVariable[!is.na(envVariable[])]+1
\tif (resistance == TRUE)
\t\t{
\t\t\tenvVariable[!is.na(envVariable[])] = 1/envVariable[!is.na(envVariable[])]
\t\t}
\tvMin = min(envVariable[!is.na(envVariable[])])
\tvMax = max(envVariable[!is.na(envVariable[])])\t\t
\t# options(digits.secs=8)
\tbuffer = ancestPosition
\tancestPosition = matrix(nrow=1, ncol=2)
\tancestPosition[1,1] = buffer[1]
\tancestPosition[1,2] = buffer[2]
\tif (extractionOfValuesOnMatrix == TRUE)
\t\t{
\t\t\tmat = raster::as.matrix(envVariable)
\t\t\txMin = extent(envVariable)@xmin
\t\t\tyMin = extent(envVariable)@ymin
\t\t\txMax = extent(envVariable)@xmax
\t\t\tyMax = extent(envVariable)@ymax
\t\t\tancestPosition[1,1] = (ancestPosition[1,1]-xMin)/(xMax-xMin)
\t\t\tancestPosition[1,1] = ceiling((ancestPosition[1,1]*(dim(mat)[2])))
\t\t\tancestPosition[1,2] = (ancestPosition[1,2]-yMin)/(yMax-yMin)
\t\t\tancestPosition[1,2] = ceiling((ancestPosition[1,2]*(dim(mat)[1])))
\t\t\textractValueOnMatrix = function(coords)
\t\t\t\t{
\t\t\t\t\tcoords_ceiling = coords
\t\t\t\t\tcoords_ceiling[1,1] = ceiling(coords_ceiling[1,1])
\t\t\t\t\tcoords_ceiling[1,2] = ceiling(dim(mat)[1]-coords_ceiling[1,2])
\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\tif (coords_ceiling[1,1] <= 0) onTheGrid = FALSE
\t\t\t\t\tif (coords_ceiling[1,2] <= 0) onTheGrid = FALSE
\t\t\t\t\tif (coords_ceiling[1,1] > dim(mat)[2]) onTheGrid = FALSE
\t\t\t\t\tif (coords_ceiling[1,2] > dim(mat)[1]) onTheGrid = FALSE
\t\t\t\t\tif (onTheGrid == FALSE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tv = NA
\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\tv = mat[coords_ceiling[1,2],coords_ceiling[1,1]]
\t\t\t\t\t\t}
\t\t\t\t\treturn(v)
\t\t\t\t}
\t\t}\telse\t{\t
\t\t\textractValueOnRaster = function(coords)
\t\t\t\t{
\t\t\t\t\treturn(raster::extract(envVariable, coords))
\t\t\t\t}\t
\t\t}
\tparticules = list()
\tnewParticule = cbind(ancestPosition,1,0,1,0)
\tparticules[[length(particules)+1]] = newParticule
\tif (showingPlots == TRUE)
\t\t{
\t\t\thisto = hist(envVariable0[!is.na(envVariable0[])], plot=F)
\t\t\tif (sum(histo$counts>0) == 2)
\t\t\t\t{
\t\t\t\t\tcols1 = colorRampPalette(brewer.pal(9,"YlOrBr"))(9)[c(1,7)]
\t\t\t\t}\telse\t{
\t\t\t\t\tcols1 = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
\t\t\t\t}
\t\t\t# if (extractionOfValuesOnMatrix == FALSE)
\t\t\t\t# {
\t\t\t\t\t# plotRaster(envVariable0, cols=cols1)
\t\t\t\t# }\telse\t{
\t\t\t\t\t# rast = raster(mat0,xmn=0,xmx=dim(mat0)[2],ymn=0,ymx=dim(mat0)[1])
\t\t\t\t\t# plotRaster(rast, cols=cols1, legend=T)
\t\t\t\t# }
\t\t\t# points(cbind(particules[[1]][1],particules[[1]][2]), pch=16, cex=0.8, col="blue")
\t\t}
\tT = startingYear + timeIntervale
\tfor (t in seq(startingYear,samplingWindow[2],timeSlice))
\t\t{
\t\t\tstatut = cbind(length(particules), t)
\t\t\tcolnames(statut) = c(); cat(statut); cat("\\n")
\t\t\tL = length(particules)
\t\t\tfor (i in 1:L)
\t\t\t\t{
\t  \t\t\t\tif (particules[[i]][5] == 1)
\t  \t\t\t\t\t{
\t   \t\t\t\t\t\tr = runif(1,0,1)
\t   \t \t\t\t\t\tif ((r < birthRate*timeSlice)|(L == 1))
\t   \t \t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tnewParticule1 = cbind(particules[[i]][1],particules[[i]][2],length(particules)+1,particules[[i]][3],1,0)
\t\t\t\t\t\t\t\t\tnewParticule2 = cbind(particules[[i]][1],particules[[i]][2],length(particules)+2,particules[[i]][3],1,0)
\t\t\t\t\t\t\t\t\tparticules[[length(particules)+1]] = newParticule1; particules[[length(particules)+1]] = newParticule2
\t\t\t\t\t\t\t\t\tparticules[[i]][5] = 0; particules[[i]][6] = t
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tfor (i in 1:length(particules))
\t\t\t\t{
\t\t\t\t\tif (particules[[i]][5] == 1)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\twhile(onTheGrid == FALSE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tcoords = cbind(particules[[i]][1],particules[[i]][2])
\t\t\t\t\t\t\t\t\tif (extractionOfValuesOnMatrix == TRUE) v1 = extractValueOnMatrix(coords)
\t\t\t\t\t\t\t\t\tif (extractionOfValuesOnMatrix == FALSE) v1 = extractValueOnRaster(coords)
\t\t\t\t\t\t\t\t\tv2 = v1/vMax
\t\t\t\t\t\t\t\t\tif (extractionOfValuesOnMatrix == TRUE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tsdX = v2*1*scalingValue
\t\t\t\t\t\t\t\t\t\t\tsdY = v2*1*scalingValue
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tsdX = v2*xres(envVariable)*scalingValue
\t\t\t\t\t\t\t\t\t\t\tsdY = v2*yres(envVariable)*scalingValue
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tdX = rnorm(1,0,sdX)
\t\t\t\t\t\t\t\t\tdY = rnorm(1,0,sdY)
\t\t\t\t\t\t\t\t\tcoords = cbind(particules[[i]][1]+dX,particules[[i]][2]+dY)
\t\t\t\t\t\t\t\t\tif (extractionOfValuesOnMatrix == TRUE) v = extractValueOnMatrix(coords)
\t\t\t\t\t\t\t\t\tif (extractionOfValuesOnMatrix == FALSE) v = extractValueOnRaster(coords)
\t\t\t\t\t\t\t\t\tif (is.na(v) == FALSE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tparticules[[i]][1] = particules[[i]][1]+dX
\t\t\t\t\t\t\t\t\t\t\tparticules[[i]][2] = particules[[i]][2]+dY
\t\t\t\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tif (t>samplingWindow[1])
\t\t\t\t{
\t\t\t\t\tfor (i in 1:length(particules))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (particules[[i]][5] == 1)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tr = runif(1,0,1)
\t\t\t\t\t\t\t\t\tif (r < (samplingRate*timeSlice))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tparticules[[i]][5] = 0; particules[[i]][6] = t
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tif (t >= T)
\t\t\t\t{
\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\t# if (extractionOfValuesOnMatrix == FALSE)
\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t# plotRaster(envVariable0, cols=cols1, new=F)
\t\t\t\t\t\t\t\t# }\telse\t{
\t\t\t\t\t\t\t\t\t# rast = raster(mat0,xmn=0,xmx=dim(mat0)[2],ymn=0,ymx=dim(mat0)[1])
\t\t\t\t\t\t\t\t\t# plotRaster(rast, cols=cols1, new=F)
\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t# coords1 = c(); coords2 = c()
\t\t\t\t\t\t\t# for (i in 1:length(particules))
\t\t\t\t\t\t\t\t# {
\t\t\t\t\t\t\t\t\t# if (particules[[i]][5] == 0) coords1 = rbind(coords1, cbind(particules[[i]][1],particules[[i]][2]))
\t\t\t\t\t\t\t\t\t# if (particules[[i]][5] == 1) coords2 = rbind(coords2, cbind(particules[[i]][1],particules[[i]][2]))
\t\t\t\t\t\t\t\t# }
\t\t\t\t\t\t\t# points(coords1, pch=16, cex=0.5, col="green3")
\t\t\t\t\t\t\t# points(coords2, pch=16, cex=0.5, col="black")
\t\t\t\t\t\t}
\t\t\t\t\tT = t + timeIntervale
\t\t\t\t}
\t\t}
\tif (extractionOfValuesOnMatrix == TRUE)
\t\t{
\t\t\tfor (i in 1:length(particules))
\t\t\t\t{
\t\t\t\t\tparticules[[i]][1] = ((particules[[i]][1]/dim(mat)[2])*(xMax-xMin))+xMin
\t\t\t\t\tparticules[[i]][2] = ((particules[[i]][2]/dim(mat)[1])*(yMax-yMin))+yMin
\t\t\t\t}
\t\t}\t
\tcsv = c()
\tfor (i in 2:length(particules))
\t\t{
\t\t\tnode1 = particules[[i]][4]; node2 = i
\t\t\tstartLon = particules[[node1]][1]; startLat = particules[[node1]][2]
\t\t\tendLon = particules[[node2]][1]; endLat = particules[[node2]][2]
\t\t\tstartYear = particules[[node1]][6]; endYear = particules[[node2]][6]
\t\t\tif ((particules[[node2]][6] == 0) & (particules[[node2]][5] == 1)) endYear = samplingWindow[2]
\t\t\tstartNodeL = startYear-startingYear; endNodeL = endYear-startingYear
\t\t\tx1 = cbind(startLon,startLat); x2 = cbind(endLon,endLat)
\t\t\tlength = endYear-startYear
\t\t\tgreatCircleDist_km = rdist.earth(x1, x2, miles=FALSE, R=NULL)
\t\t\ttreeID = length; treeID[] = -9999
\t\t\tcsv = rbind(csv, cbind(node1,node2,startLat,startLon,endLat,endLon,endNodeL,startNodeL,startYear,endYear,greatCircleDist_km,treeID))
\t\t}
\tcolNames = c("node1","node2","startLat","startLon","endLat","endLon","endNodeL","startNodeL",
\t\t\t\t\t"startYear","endYear","greatCircleDist_km","treeID")
\tcolnames(csv) = colNames
\tcsv = as.matrix(csv); row.names(csv) = c()
\tunsampledTipNodes = as.numeric(which(csv[,"endYear"] == samplingWindow[2]))
\tif (length(unsampledTipNodes) > 1)
\t\t{
\t\t\tunsampledTipNodes = sample(unsampledTipNodes, length(unsampledTipNodes)-1, replace=F)
\t\t\tindices = c(1:dim(csv)[1]); indices = indices[!indices%in%unsampledTipNodes]
\t\t\tcsv = csv[indices,]
\t\t\tbranchesToRemove = c()
\t\t\tfor (j in 1:dim(csv)[1])
\t\t\t\t{
\t\t\t\t\tif (length(which(csv[,"node1"]==csv[j,"node2"])) == 1)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tindex = which(csv[,"node1"]==csv[j,"node2"])
\t\t\t\t\t\t\tcsv[j,"node2"] = csv[index,"node2"]
\t\t\t\t\t\t\tcsv[j,"endLat"] = csv[index,"endLat"]
\t\t\t\t\t\t\tcsv[j,"endNodeL"] = csv[index,"endNodeL"]
\t\t\t\t\t\t\tcsv[j,"endYear"] = csv[index,"endYear"]
\t\t\t\t\t\t\t# x1 = cbind(csv[j,"startLon"],csv[j,"startLat"]); x2 = cbind(csv[j,"endLon"],csv[j,"endLat"])
\t\t\t\t\t\t\t# csv[j,"greatCircleDist_km"] = rdist.earth(x1, x2, miles=FALSE, R=NULL)
\t\t\t\t\t\t\tcsv[j,"greatCircleDist_km"] = NaN
\t\t\t\t\t\t\tbranchesToRemove = c(branchesToRemove, index)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tindices = c(1:dim(csv)[1]); indices = indices[!indices%in%branchesToRemove]
\t\t\tcsv = csv[indices,]
\t\t}
\tif (dim(csv)[1] > 2)
\t\t{
\t\t\tisolatedClade = TRUE
\t\t\twhile (isolatedClade == TRUE)
\t\t\t\t{
\t\t\t\t\tbranchesToRemove = c()
\t\t\t\t\tfor (j in 3:dim(csv)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (!csv[j,"node1"]%in%csv[,"node2"])
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tbranchesToRemove = c(branchesToRemove, j)
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (length(branchesToRemove) == 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tisolatedClade = FALSE
\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\tindices = c(1:dim(csv)[1]); indices = indices[!indices%in%branchesToRemove]
\t\t\t\t\t\t\tcsv = csv[indices,]\t
\t\t\t\t\t\t}
\t\t\t\t}
\t\t}
\tbranches1 = c(); nodes1 = c()
\ttipNodes = c()
\tfor (j in 1:dim(csv)[1])
\t\t{
\t\t\ttipNode = TRUE
\t\t\tfor (k in 1:dim(csv)[1])
\t\t\t\t{
\t\t\t\t\tif (csv[j,"node2"] == csv[k,"node1"]) tipNode = FALSE
\t\t\t\t}
\t\t\tif (tipNode == TRUE)
\t\t\t\t{
\t\t\t\t\tlength = csv[j,"endYear"]-csv[j,"startYear"]
\t\t\t\t\tbranches1 = c(branches1, paste(csv[j,"node2"],":",length,sep=""))
\t\t\t\t\tnodes1 =  c(nodes1, csv[j,"node2"])
\t\t\t\t\ttipNodes = rbind(tipNodes, c(csv[j,"node2"], csv[j,"startYear"]))
\t\t\t\t}
\t\t}
\tcoalescenceEvents = length(branches1)-1
\tc = 0; 
\twhile (length(branches1) != 1)
\t\t{
\t\t\tc = c+1
\t\t\tbranches2 = branches1; nodes2 = nodes1
\t\t\tnodesToRemove = c()
\t\t\tfor (j in 2:length(branches1))
\t\t\t\t{
\t\t\t\t\tnodeA = nodes1[j]
\t\t\t\t\tfor (k in 1:(j-1))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tnodeB = nodes1[k]
\t\t\t\t\t\t\tif (csv[(csv[,"node2"]==nodeA),"node1"] == csv[(csv[,"node2"]==nodeB),"node1"])
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tcoalescenceEvents = coalescenceEvents-1
\t\t\t\t\t\t\t\t\tnodesToRemove = c(nodesToRemove, k, j)
\t\t\t\t\t\t\t\t\tnode2 = csv[(csv[,"node2"]==csv[(csv[,"node2"]==nodeA),"node1"]),"node2"]
\t\t\t\t\t\t\t\t\tlength = csv[(csv[,"node2"]==node2),"endYear"]-csv[(csv[,"node2"]==node2),"startYear"]
\t\t\t\t\t\t\t\t\tif (coalescenceEvents > 0)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tbranches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",":",length,sep=""))
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tbranches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",";",length,sep=""))
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tnodes2 = c(nodes2, node2)
\t\t\t\t\t\t\t\t\t# print(paste(coalescenceEvents,nodeB,nodeA,node2,sep=" "))
\t\t\t\t\t\t\t\t\t# print(c(node2,j,k))
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tbranches2 = branches2[-nodesToRemove]
\t\t\tnodes2 = nodes2[-nodesToRemove]
\t\t\tbranches1 = branches2; nodes1 = nodes2
\t\t}
\tbranches3 = branches1
\ttree = read.tree(text=branches1)
\tcsv = csv[order(csv[,"startYear"],decreasing=F),]
\tcsv1 = csv[1,]; csv2 = csv[2:dim(csv)[1],]
\tcsv2 = csv2[order(csv2[,"endYear"],decreasing=F),]
\tcsv = rbind(csv1, csv2)
\tif (showingPlots == TRUE)
\t\t{
\t\t\tcol_start = colorRampPalette(brewer.pal(9,'PuBu'))(101)[1]
\t\t\tcols2 = colorRampPalette(brewer.pal(9,'PuBu'))(101)[(((csv[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
\t\t\tplotRaster(envVariable0, cols=cols1, addLegend=TRUE)
\t\t\tfor (i in 1:dim(csv)[1])
\t\t\t\t{
\t\t\t\t\tsegments(csv[i,"startLon"], csv[i,"startLat"], csv[i,"endLon"], csv[i,"endLat"], col=rgb(1,1,1,255,maxColorValue=255), lwd=0.5)
\t\t\t\t}
\t\t\tfor (i in 1:dim(csv)[1])
\t\t\t\t{
\t\t\t\t\tif (i == 1)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tpoints(cbind(csv[i,"startLon"],csv[i,"startLat"]), pch=16, cex=0.9, col=col_start)
\t\t\t\t\t\t\tpoints(cbind(csv[i,"startLon"],csv[i,"startLat"]), pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
\t\t\t\t\t\t}
\t\t\t\t\tpoints(cbind(csv[i,"endLon"],csv[i,"endLat"]), pch=16, cex=0.9, col=cols2[i])
\t\t\t\t\tpoints(cbind(csv[i,"endLon"],csv[i,"endLat"]), pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
\t\t\t\t}
\t\t\tdev.new()
\t\t\tcol_start = colorRampPalette(brewer.pal(9,'PuBu'))(101)[1]
\t\t\tcols3 = colorRampPalette(brewer.pal(9,'PuBu'))(101)[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
\t\t\tplot(tree, show.tip.label=F, edge.width=0.5)
\t\t\tfor (i in 1:dim(tree$edge)[1])
\t\t\t\t{
\t\t\t\t\tif (i == 1)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tnodelabels(node=tree$edge[i,1], pch=16, cex=0.9, col=col_start)
\t\t\t\t\t\t\tnodelabels(node=tree$edge[i,1], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
\t\t\t\t\t\t}
\t\t\t\t\tnodelabels(node=tree$edge[i,2], pch=16, cex=0.9, col=cols3[i])
\t\t\t\t\tnodelabels(node=tree$edge[i,2], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
\t\t\t\t}
\t\t}
\toutputs = list(csv, tree)
\treturn(outputs)
}
