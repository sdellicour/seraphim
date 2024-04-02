simulatorRRW3 = function(envVariable=raster(matrix(c(runif(600,5,10),runif(1000,0,5)),nrow=40,ncol=40)),
						 resistance=TRUE, scalingValue=1, ancestPosition=c(0.4,0.5), birthRate=0.2,
						 samplingRate=0.2, startingYear=0, samplingWindow=c(10,50), timeSlice=0.1,
						 timeIntervale=1, showingPlots=FALSE, extractionOfValuesOnMatrix=FALSE) {
	
	envVariable0 = envVariable
	mat0 = raster::as.matrix(envVariable)
	envVariable[!is.na(envVariable[])] = envVariable[!is.na(envVariable[])]+1
	if (resistance == TRUE)
		{
			envVariable[!is.na(envVariable[])] = 1/envVariable[!is.na(envVariable[])]
		}
	vMin = min(envVariable[!is.na(envVariable[])])
	vMax = max(envVariable[!is.na(envVariable[])])		
	# options(digits.secs=8)
	buffer = ancestPosition
	ancestPosition = matrix(nrow=1, ncol=2)
	ancestPosition[1,1] = buffer[1]
	ancestPosition[1,2] = buffer[2]
	if (extractionOfValuesOnMatrix == TRUE)
		{
			mat = raster::as.matrix(envVariable)
			xMin = extent(envVariable)@xmin
			yMin = extent(envVariable)@ymin
			xMax = extent(envVariable)@xmax
			yMax = extent(envVariable)@ymax
			ancestPosition[1,1] = (ancestPosition[1,1]-xMin)/(xMax-xMin)
			ancestPosition[1,1] = ceiling((ancestPosition[1,1]*(dim(mat)[2])))
			ancestPosition[1,2] = (ancestPosition[1,2]-yMin)/(yMax-yMin)
			ancestPosition[1,2] = ceiling((ancestPosition[1,2]*(dim(mat)[1])))
			extractValueOnMatrix = function(coords)
				{
					coords_ceiling = coords
					coords_ceiling[1,1] = ceiling(coords_ceiling[1,1])
					coords_ceiling[1,2] = ceiling(dim(mat)[1]-coords_ceiling[1,2])
					onTheGrid = TRUE
					if (coords_ceiling[1,1] <= 0) onTheGrid = FALSE
					if (coords_ceiling[1,2] <= 0) onTheGrid = FALSE
					if (coords_ceiling[1,1] > dim(mat)[2]) onTheGrid = FALSE
					if (coords_ceiling[1,2] > dim(mat)[1]) onTheGrid = FALSE
					if (onTheGrid == FALSE)
						{
							v = NA
						}	else	{
							v = mat[coords_ceiling[1,2],coords_ceiling[1,1]]
						}
					return(v)
				}
		}	else	{	
			extractValueOnRaster = function(coords)
				{
					return(raster::extract(envVariable, coords))
				}	
		}
	particules = list()
	newParticule = cbind(ancestPosition,1,0,1,0)
	particules[[length(particules)+1]] = newParticule
	if (showingPlots == TRUE)
		{
			histo = hist(envVariable0[!is.na(envVariable0[])], plot=F)
			if (sum(histo$counts>0) == 2)
				{
					cols1 = colorRampPalette(brewer.pal(9,"YlOrBr"))(9)[c(1,7)]
				}	else	{
					cols1 = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
				}
		}
	T = startingYear + timeIntervale
	for (t in seq(startingYear,samplingWindow[2],timeSlice))
		{
			statut = cbind(length(particules), t)
			colnames(statut) = c(); cat(statut); cat("\n")
			L = length(particules)
			for (i in 1:L)
				{
	  				if (particules[[i]][5] == 1)
	  					{
	   						r = runif(1,0,1)
	   	 					if ((r < birthRate*timeSlice)|(L == 1))
	   	 						{
									newParticule1 = cbind(particules[[i]][1],particules[[i]][2],length(particules)+1,particules[[i]][3],1,0)
									newParticule2 = cbind(particules[[i]][1],particules[[i]][2],length(particules)+2,particules[[i]][3],1,0)
									particules[[length(particules)+1]] = newParticule1; particules[[length(particules)+1]] = newParticule2
									particules[[i]][5] = 0; particules[[i]][6] = t
								}
						}
				}
			for (i in 1:length(particules))
				{
					if (particules[[i]][5] == 1)
						{
							onTheGrid = FALSE; c = 0
							while (onTheGrid == FALSE)
								{
									c = c+1
									if (c == 1000) print(i)
									coords = cbind(particules[[i]][1],particules[[i]][2])
									if (extractionOfValuesOnMatrix == TRUE) v1 = extractValueOnMatrix(coords)
									if (extractionOfValuesOnMatrix == FALSE) v1 = extractValueOnRaster(coords)
									v2 = v1/vMax
									if (extractionOfValuesOnMatrix == TRUE)
										{
											sdX = v2*1*scalingValue
											sdY = v2*1*scalingValue
										}	else	{
											sdX = v2*xres(envVariable)*scalingValue
											sdY = v2*yres(envVariable)*scalingValue
										}
									dX = rnorm(1,0,sdX)
									dY = rnorm(1,0,sdY)
									coords = cbind(particules[[i]][1]+dX,particules[[i]][2]+dY)
									if (extractionOfValuesOnMatrix == TRUE) v = extractValueOnMatrix(coords)
									if (extractionOfValuesOnMatrix == FALSE) v = extractValueOnRaster(coords)
									if (is.na(v) == FALSE)
										{
											particules[[i]][1] = particules[[i]][1]+dX
											particules[[i]][2] = particules[[i]][2]+dY
											onTheGrid = TRUE
										}
								}
						}
				}
			if (t>samplingWindow[1])
				{
					for (i in 1:length(particules))
						{
							if (particules[[i]][5] == 1)
								{
									r = runif(1,0,1)
									if (r < (samplingRate*timeSlice))
										{
											particules[[i]][5] = 0; particules[[i]][6] = t
										}
								}
						}
				}
			if (t >= T)
				{
					T = t + timeIntervale
				}
		}
	if (extractionOfValuesOnMatrix == TRUE)
		{
			for (i in 1:length(particules))
				{
					particules[[i]][1] = ((particules[[i]][1]/dim(mat)[2])*(xMax-xMin))+xMin
					particules[[i]][2] = ((particules[[i]][2]/dim(mat)[1])*(yMax-yMin))+yMin
				}
		}	
	csv = c()
	for (i in 2:length(particules))
		{
			node1 = particules[[i]][4]; node2 = i
			startLon = particules[[node1]][1]; startLat = particules[[node1]][2]
			endLon = particules[[node2]][1]; endLat = particules[[node2]][2]
			startYear = particules[[node1]][6]; endYear = particules[[node2]][6]
			if ((particules[[node2]][6] == 0) & (particules[[node2]][5] == 1)) endYear = samplingWindow[2]
			startNodeL = startYear-startingYear; endNodeL = endYear-startingYear
			x1 = cbind(startLon,startLat); x2 = cbind(endLon,endLat)
			length = endYear-startYear
			greatCircleDist_km = rdist.earth(x1, x2, miles=FALSE, R=NULL)
			treeID = length; treeID[] = -9999
			csv = rbind(csv, cbind(node1,node2,startLat,startLon,endLat,endLon,endNodeL,startNodeL,startYear,endYear,greatCircleDist_km,treeID))
		}
	colNames = c("node1","node2","startLat","startLon","endLat","endLon","endNodeL","startNodeL","startYear","endYear","greatCircleDist_km","treeID")
	colnames(csv) = colNames
	csv = as.matrix(csv); row.names(csv) = c()
	unsampledTipNodes = as.numeric(which(csv[,"endYear"] == samplingWindow[2]))
	if (length(unsampledTipNodes) > 1)
		{
			unsampledTipNodes = sample(unsampledTipNodes, length(unsampledTipNodes)-1, replace=F)
			indices = c(1:dim(csv)[1]); indices = indices[!indices%in%unsampledTipNodes]
			csv = csv[indices,]
			branchesToRemove = c()
			for (j in 1:dim(csv)[1])
				{
					if (length(which(csv[,"node1"]==csv[j,"node2"])) == 1)
						{
							index = which(csv[,"node1"]==csv[j,"node2"])
							csv[j,"node2"] = csv[index,"node2"]
							csv[j,"endLat"] = csv[index,"endLat"]
							csv[j,"endNodeL"] = csv[index,"endNodeL"]
							csv[j,"endYear"] = csv[index,"endYear"]
							# x1 = cbind(csv[j,"startLon"],csv[j,"startLat"]); x2 = cbind(csv[j,"endLon"],csv[j,"endLat"])
							# csv[j,"greatCircleDist_km"] = rdist.earth(x1, x2, miles=FALSE, R=NULL)
							csv[j,"greatCircleDist_km"] = NaN
							branchesToRemove = c(branchesToRemove, index)
						}
				}
			indices = c(1:dim(csv)[1]); indices = indices[!indices%in%branchesToRemove]
			csv = csv[indices,]
		}
	if (dim(csv)[1] > 2)
		{
			isolatedClade = TRUE
			while (isolatedClade == TRUE)
				{
					branchesToRemove = c()
					for (j in 3:dim(csv)[1])
						{
							if (!csv[j,"node1"]%in%csv[,"node2"])
								{
									branchesToRemove = c(branchesToRemove, j)
								}
						}
					if (length(branchesToRemove) == 0)
						{
							isolatedClade = FALSE
						}	else	{
							indices = c(1:dim(csv)[1]); indices = indices[!indices%in%branchesToRemove]
							csv = csv[indices,]	
						}
				}
		}
	branches1 = c(); nodes1 = c()
	tipNodes = c()
	for (j in 1:dim(csv)[1])
		{
			tipNode = TRUE
			for (k in 1:dim(csv)[1])
				{
					if (csv[j,"node2"] == csv[k,"node1"]) tipNode = FALSE
				}
			if (tipNode == TRUE)
				{
					length = csv[j,"endYear"]-csv[j,"startYear"]
					branches1 = c(branches1, paste(csv[j,"node2"],":",length,sep=""))
					nodes1 =  c(nodes1, csv[j,"node2"])
					tipNodes = rbind(tipNodes, c(csv[j,"node2"], csv[j,"startYear"]))
				}
		}
	coalescenceEvents = length(branches1)-1
	c = 0; 
	while (length(branches1) != 1)
		{
			c = c+1
			branches2 = branches1; nodes2 = nodes1
			nodesToRemove = c()
			for (j in 2:length(branches1))
				{
					nodeA = nodes1[j]
					for (k in 1:(j-1))
						{
							nodeB = nodes1[k]
							if (csv[(csv[,"node2"]==nodeA),"node1"] == csv[(csv[,"node2"]==nodeB),"node1"])
								{
									coalescenceEvents = coalescenceEvents-1
									nodesToRemove = c(nodesToRemove, k, j)
									node2 = csv[(csv[,"node2"]==csv[(csv[,"node2"]==nodeA),"node1"]),"node2"]
									length = csv[(csv[,"node2"]==node2),"endYear"]-csv[(csv[,"node2"]==node2),"startYear"]
									if (coalescenceEvents > 0)
										{
											branches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",":",length,sep=""))
										}	else	{
											branches2 = c(branches2, paste("(",branches1[k],",",branches1[j],")",";",length,sep=""))
										}
									nodes2 = c(nodes2, node2)
								}
						}
				}
			branches2 = branches2[-nodesToRemove]
			nodes2 = nodes2[-nodesToRemove]
			branches1 = branches2; nodes1 = nodes2
		}
	branches3 = branches1
	tree = read.tree(text=branches1)
	csv = csv[order(csv[,"startYear"],decreasing=F),]
	csv1 = csv[1,]; csv2 = csv[2:dim(csv)[1],]
	csv2 = csv2[order(csv2[,"endYear"],decreasing=F),]
	csv = rbind(csv1, csv2)
	if (showingPlots == TRUE)
		{
			col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
			cols2 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((csv[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
			plotRaster(envVariable0, cols=cols1, addLegend=TRUE)
			for (i in 1:dim(csv)[1])
				{
					segments(csv[i,"startLon"], csv[i,"startLat"], csv[i,"endLon"], csv[i,"endLat"], col=rgb(1,1,1,255,maxColorValue=255), lwd=0.5)
				}
			for (i in 1:dim(csv)[1])
				{
					if (i == 1)
						{
							points(cbind(csv[i,"startLon"],csv[i,"startLat"]), pch=16, cex=0.9, col=col_start)
							points(cbind(csv[i,"startLon"],csv[i,"startLat"]), pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
						}
					points(cbind(csv[i,"endLon"],csv[i,"endLat"]), pch=16, cex=0.9, col=cols2[i])
					points(cbind(csv[i,"endLon"],csv[i,"endLat"]), pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
				}
			dev.new()
			col_start = colorRampPalette(brewer.pal(9,"PuBu"))(101)[1]
			cols3 = colorRampPalette(brewer.pal(9,"PuBu"))(101)[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
			plot(tree, show.tip.label=F, edge.width=0.5)
			for (i in 1:dim(tree$edge)[1])
				{
					if (i == 1)
						{
							nodelabels(node=tree$edge[i,1], pch=16, cex=0.9, col=col_start)
							nodelabels(node=tree$edge[i,1], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
						}
					nodelabels(node=tree$edge[i,2], pch=16, cex=0.9, col=cols3[i])
					nodelabels(node=tree$edge[i,2], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
				}
		}
	outputs = list(csv, tree)
	return(outputs)
}
