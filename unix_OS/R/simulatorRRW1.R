simulatorRRW1 <-
function(tree, rates, sigmas=c(0.1,0.1), cor=0, envVariables=list(), mostRecentSamplingDatum,
						 ancestPosition, reciprocalRates=TRUE, n1=100, n2=0, showingPlots=FALSE, newPlot=TRUE, fixedNodes=c()) {
	
	rotation = function(pt1, pt2, angle)
		{
			s = sin(angle); c = cos(angle)
			x = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
			x_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
			x_new = x_new+pt1[1]; y_new = y_new+pt1[2]
			return(c(x_new,y_new))
		}
	nodesOnly = FALSE; pointCol = "red"
	colNames= c("node1","node2","length","startLat","startLon","endLat","endLon",
				"endNodeL","startNodeL","startYear","endYear","greatCircleDist_km")
	simulation = matrix(nrow=length(tree$edge.length), ncol=length(colNames))
	colnames(simulation) = colNames
	# if (model == "fixed") phi_b = rep(1, length(tree$edge.length))
	# if (model == "cauchy") phi_b = rgamma(length(tree$edge.length), shape=0.5, scale=0.5)
	# if (model == "gamma") phi_b = rgamma(length(tree$edge.length), shape=halfDF, scale=halfDF)
	# if (model == "logN") phi_b = rlnorm(length(tree$edge.length), meanlog=1, sdlog=sdLogN)
	# if (model == "cauchy") sd_BM = sqrt(tree$edge.length/phi_b) # corresponds to BEAST reciprocalRates="true"
	# if (model == "gamma") sd_BM = sqrt(tree$edge.length/phi_b) # corresponds to BEAST reciprocalRates="true"
	# if (model == "logN") sd_BM = sqrt(tree$edge.length*phi_b) # corresponds to BEAST reciprocalRates="false"
	phi_b = rates
	if (reciprocalRates == FALSE) sd_BM = sqrt(tree$edge.length*phi_b)
	if (reciprocalRates == TRUE) sd_BM = sqrt(tree$edge.length/phi_b)
	nd = node.depth(tree); nd_max = max(nd)
	t1 = rep(ancestPosition[2], length(tree$tip.label)+tree$Nnode)
	t2 = rep(ancestPosition[1], length(tree$tip.label)+tree$Nnode)
	if (showingPlots == TRUE)
		{
			if ((newPlot == TRUE)&(length(envVariables) > 0))
				{
					plotRaster(rast=envVariables[[1]], cols="gray90", colNA="white", addBox=F)
				}
			points(cbind(ancestPosition[1],ancestPosition[2]), pch=16, col=pointCol, cex=0.5)
		}
	i = 0
	while (i != (nd_max-1))
		{
			i = i+1
			my_nodes = which(nd==nd_max-i)
			if (length(my_nodes) > 0)
				{
					for (j in 1:length(my_nodes))
						{
							my_node = my_nodes[j]; # print(c(i,my_nodes[j]))
							parent_branch = match(my_node, tree$edge[,2])
							parent_node = tree$edge[parent_branch,1] 
							simulatingNode = TRUE
							if (my_node%in%fixedNodes)
								{
									simulatingNode = FALSE
									index = which(tree$edge[,2] == my_node)
									new_t1 = tree$annotations[[index]]$location[[1]]
									new_t2 = tree$annotations[[index]]$location[[2]]
								}
							if (simulatingNode == TRUE)
								{
									onTheArea = FALSE
									sd_bm = sd_BM[parent_branch]
									increment1 = rnorm(1)*sd_bm
									increment2 = rnorm(1)*sd_bm
									increment2 = sigmas[2]*((cor*increment1)+(sqrt(1-cor^2)*increment2))
									increment1 = sigmas[1]*increment1
									new_t1 = t1[parent_node] + increment1 
									new_t2 = t2[parent_node] + increment2
									if (length(envVariables) > 0)
										{
											onTheArea = TRUE
											for (k in 1:length(envVariables))
												{
													if (is.na(raster::extract(envVariables[[k]],cbind(new_t2,new_t1))))
														{
															onTheArea = FALSE
														}
												}
										}	else	{
											onTheArea = TRUE
										}
									if (onTheArea == FALSE)
										{
											c2 = 0; c1 = 0
											pt1 = cbind(t2[parent_node],t1[parent_node])
											pt2 = cbind(new_t2,new_t1)
											while (onTheArea == FALSE)
												{
													c2 = c2+1; # print(c(c2,c1))
													if (n1 > 0)
														{
															angle = (2*pi)*runif(1)
															pt2_rotated = rotation(pt1, pt2, angle)
															onTheArea = TRUE
															for (k in 1:length(envVariables))
																{
																	if (is.na(raster::extract(envVariables[[k]],cbind(pt2_rotated[1],pt2_rotated[2]))))
																		{
																			onTheArea = FALSE
																		}	else		{
																			new_t1 = pt2_rotated[2]
																			new_t2 = pt2_rotated[1]
																		}
																}
														}
													if (c2 > n1)
														{
															c2 = 0; c1 = c1+1
															if (c1 > n2)
																{
																	onTheArea = TRUE; i = 0
																}	else	{
																	# print(paste0("...re-simulating a branch - node:", my_node))
																	increment1 = rnorm(1)*sd_bm
																	increment2 = rnorm(1)*sd_bm
																	increment2 = sigmas[2]*((cor*increment1)+(sqrt(1-cor^2)*increment2))
																	increment1 = sigmas[1]*increment1
																	new_t1 = t1[parent_node] + increment1 
																	new_t2 = t2[parent_node] + increment2
																	onTheArea = TRUE
																	for (k in 1:length(envVariables))
																		{
																			if (is.na(raster::extract(envVariables[[k]],cbind(new_t2,new_t1))))
																				{
																					onTheArea = FALSE
																				}
																		}
																}
														}
												}
										}
								}
							t1[my_node] = new_t1	
							t2[my_node] = new_t2
							if (showingPlots == TRUE)
								{
									if (nodesOnly == FALSE)
										{
											segments(t2[parent_node], t1[parent_node], new_t2, new_t1, col=pointCol, lwd=0.2)
											points(cbind(new_t2,new_t1), pch=16, col=pointCol, cex=0.25)
										}	else	{
											points(cbind(new_t2,new_t1), pch=16, col=pointCol, cex=0.25)
										}
								}
						}
					if (i == 0)
						{
							cat(paste0("...re-starting the simulation\\n"))
							t1 = rep(ancestPosition[2], length(tree$tip.label)+tree$Nnode)
							t2 = rep(ancestPosition[1], length(tree$tip.label)+tree$Nnode)
							if (showingPlots == TRUE)
								{
									plotRaster(rast=envVariables[[1]], cols="gray90", colNA="white", addBox=F, new=F)
									points(cbind(ancestPosition[1],ancestPosition[2]), pch=16, col=pointCol, cex=0.5)
								}
						}
				}
		}
	x = t2; y = t1
	for (i in 1:dim(tree$edge)[1])
		{
			node_i = tree$edge[i,1]
			node_f = tree$edge[i,2]
			simulation[i,"node1"] = node_i
			simulation[i,"node2"] = node_f
			simulation[i,"length"] = tree$edge.length[i]
			simulation[i,"startLat"] = y[node_i]
			simulation[i,"startLon"] = x[node_i]
			simulation[i,"endLat"] = y[node_f]
			simulation[i,"endLon"] = x[node_f]
			x1 = cbind(x[node_i],y[node_i]); x2 = cbind(x[node_f],y[node_f])
			simulation[i,"greatCircleDist_km"] = rdist.earth(x1, x2, miles=FALSE, R=NULL)
		}	
	l = length(simulation[,1])
	ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
	for (j in 1:l)
		{
			subMat = simulation[j,2]
			subMat = subset(simulation,simulation[,2]==subMat)
			ll[j,1] = subMat[,3]
			subMat = subMat[1,1]
			subMat1 = subset(simulation,simulation[,2]==subMat)
			for (k in 1:l)
				{
	 				if (nrow(subMat1) > 0)
	 					{
							ll[j,k+1] = subMat1[,3]
 							subMat2 = subMat1[1,1]
 							subMat1 = subset(simulation,simulation[,2]==subMat2)
 						}
 				}
		}
	endNodeL = rowSums(ll) # root to node distance for each node
	simulation[,"endNodeL"] = endNodeL
	startNodeL = matrix(1:l,nrow=l,ncol=1)
	startNodeL[] = 0
	for (j in 1:l)
		{
			r = simulation[j,1]
			s = subset(simulation,simulation[,2]==r)
			for (k in 1:l)
				{
					if (nrow(s) > 0)
						{
							startNodeL[j,1] = s[,"endNodeL"]
 						}
 				}	
	 	}
	simulation[,"startNodeL"] = startNodeL
	maxEndLIndice = which.max(simulation[,"endNodeL"])
 	maxEndL = simulation[maxEndLIndice,"endNodeL"]
 	endYear = matrix(simulation[,"endNodeL"]-maxEndL)
 	endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
	startYear = matrix(simulation[,"startNodeL"]-maxEndL)
 	startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
 	simulation[,c("startYear","endYear")] = cbind(startYear,endYear)
 	if (showingPlots == TRUE) dev.off()
	return(simulation)
}
