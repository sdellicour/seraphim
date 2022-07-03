treeExtractions = function(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName) {

	# registerDoMC(cores=nberOfCores)
	dir.create(file.path(localTreesDirectory), showWarnings=F)
	MCCtree = FALSE
	newTreeList = list(); count = 1
	for (i in 1:length(allTrees)) # to separate trees
		{	
			if (substring(allTrees[[i]],1,4) == "tree")
				{
					newTreeList[[count]] = allTrees[[i]]
					count = count+1
				}
		}
	nberOfTrees	= length(newTreeList)	
	newTreeListMinusBurnIn = newTreeList[(burnIn+1):nberOfTrees]
	regularSampledTrees = newTreeList[nberOfTreesToSample]
	if (randomSampling == TRUE)
		{
			sampledNewTreeListMinusBurnIn = sample(newTreeListMinusBurnIn,nberOfTreesToSample)
		}	else	{
			intervale = floor(length(newTreeListMinusBurnIn)/nberOfTreesToSample)
			index = 0
			for (i in 1:nberOfTreesToSample)
				{
					index = index + intervale
					regularSampledTrees[i] = newTreeListMinusBurnIn[index]
				}
			sampledNewTreeListMinusBurnIn = regularSampledTrees
		}
	
	buffer = list()
	# buffer = foreach(i = 1:length(sampledNewTreeListMinusBurnIn)) %dopar% {
	for (i in 1:length(sampledNewTreeListMinusBurnIn)) {	
			extractionsID = runif(1,0,1)
			properTreeName = paste("ProperTree_", extractionsID, ".tree", sep="")
			tree = sampledNewTreeListMinusBurnIn[i]
			if (MCCtree == FALSE)
				{
					treeID = gsub(pattern="\\[.+\\];", replacement="", tree)
					treeID = gsub(pattern=" ", replacement="", treeID)
					treeID = unlist(strsplit(treeID, "_"))[2]
					cat("Exctracting information from sampled tree n°", treeID, "\n", sep="");
				}	else		{
					treeID = "MCC_tree"
				}
			# (1) Creation of "tab", a data frame containing coordinates of each node:
			tab = gsub(pattern="tree.+\\[&R\\]", replacement="", tree)	
			tab = gsub(pattern="\\}\\]:\\[&rate=", replacement="\\},rate=", tab) 
			tab = gsub(pattern="\\}\\]:\\[&+[[:alnum:]]+[[:punct:]]+rate=", replacement="\\},rate=", tab)
			tab = gsub("\\[&CO", "\\&&CO", tab)
			tab = gsub("\\[&GTEV", "\\&&GTEV", tab)
			tab = unlist(strsplit(tab, "\\["))[-1]
			tab = gsub("&", "", tab)
			tab = gsub("\\}.+$", "", tab)		
			tab = gsub("\\{", "", tab)
			tab = gsub(paste(".*",coordinateAttributeName,"=",sep=""), "location1=", tab)
			tab = gsub(",", ",location2=", tab)
	    	foo = function(x) {x = unlist(strsplit(x, ",")); x}
	    	tab = lapply(tab, foo)
	    	tab_transit = data.frame()
	    	colNames = unique(gsub("=.+$", "", unlist(tab)))
	    	colNames = gsub(" ", "", colNames)
	    	for (j in seq(tab))
	    		{
	    			for (k in seq(colNames))
	    				{	
	    					ind = grep(paste("^",colNames[k],"=",sep=""), tab[[j]])
	    					if (length(ind) > 0)
	    						{
	    							v = as.numeric(gsub(paste(colNames[k],"=",sep=""),"",tab[[j]][ind]))
	    							tab_transit[j,k] = v
	    						}
	    				}
	    		}
	    	tab = tab_transit
	    	colnames(tab) = colNames   
	    	all = which(!is.na(tab[,1]))
	    	 
			# (2) Removing of all the spatial & mutation rate information:
			tree = gsub("\\[[^]]*\\]", "", tree) # remove all the "[...]"
			tree = gsub("tree STATE_.+[[:space:]]+=[[:space:]]+", "", tree)
			write(tree, file=properTreeName)
			
			# (3) Creation of "mat", a matrix containing coordinates of each node:
				# tree1 will contain node names (if tip nodes) and branch lengths
			tree1 = unlist(strsplit(tree, "\\("))
	   		tree1 = unlist(strsplit(tree1, ","))
	   		tree1 = unlist(strsplit(tree1, "\\)"))
		   		# tree2 will only contain all the branch lengths
			tree2 = unlist(strsplit(tree, ":"))[-1]
			tree2 = gsub("[( | ) | ;]", "", tree2)	
			tree2 = strsplit(tree2, ",")
			foo = function(x) x = head(x, 1)
			tree2 = unlist(lapply(tree2, foo))
			tree2 = paste("", tree2, sep=":")
			tree2 = c(tree2, ";")
			allStats = vector(mode="list", length=dim(tab)[2])
			for (j in seq(allStats))
				{
					tree3 = tree
					val = tab[,j]
					tree4 = paste(val, tree2, sep="")
						# tree4 will only contain one coordinate and the branch length
					for (k in all)
						{
							# tree3 = gsub(tree1[k], tree4[k], tree3) # PROBLEM
							tree3 = gsub(paste("\\(",tree1[k],sep=""), paste("\\(",tree4[k],sep=""), tree3)
							tree3 = gsub(paste("\\)",tree1[k],sep=""), paste("\\)",tree4[k],sep=""), tree3)
							tree3 = gsub(paste(",",tree1[k],sep=""), paste(",",tree4[k],sep=""), tree3)
								# node/tip names are replaced by the coordinate
	 						tree5 = read.tree(text = tree3)
							tree6 = c(tree5$tip.label, tree5$node.label)
							tree6[tree6 == "NA"] = 9999
	 						tree6 = as.numeric(tree6)
							tree6[tree6 == 9999] = NA
							allStats[[j]] = tree6
							names(allStats)[j] = colnames(tab)[j]	
						}	# allStats will contain of list of the two coordinates
				}			# but in the correct order (i.e. the "tree$edge order")
			mat = matrix(unlist(allStats), ncol=length(allStats), byrow=F)
	
			# (4) Adding the spatial information to the matrix:
			tree = read.tree(properTreeName)
			treeEdges = tree$edge
			treeEdgesLengths = tree$edge.length
			treeEdgesLengths = as.matrix(treeEdgesLengths,nrows=(b[1]),ncol=1)
			treeEdgesLengths = cbind(treeEdges,treeEdgesLengths)
			mat = cbind(mat,c(1:length(mat[,1])))  
			mat = cbind(mat,c(1:length(mat[,1]))) 
			colnames(mat) = c("latitude","longitude","node1","node2")
			colnames(treeEdgesLengths) = c("node1","node2","length")
			mat_transit = merge(treeEdgesLengths,mat,"node2")
			colnames(mat_transit) = c("node2","node1","length","endStates","endRate","x")
			mat_transit = merge(mat_transit,mat,"node1")
			colnames(mat_transit) = c("node1","node2","length","endLat","endLon","x","startLat","startLon","y")
			mat_transit$x = NULL; mat_transit$y = NULL
			mat = mat_transit[,c(1,2,3,6,7,4,5)] # changes the order of the columns
	
			# (5) Adding the temporal information to the matrix:
			l = length(mat[,1])
			if (l > 0)
				{
					ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
					for (j in 1:l)
						{
							subMat = mat[j,2]
							subMat = subset(mat,mat[,2]==subMat)
							ll[j,1] = subMat[,3]
							subMat = subMat[1,1]
							subMat1 = subset(mat,mat[,2]==subMat)
							for (k in 1:l)
								{
		 							if (nrow(subMat1) > 0)
		 								{
											ll[j,k+1] = subMat1[,3]
	 										subMat2 = subMat1[1,1]
	 										subMat1 = subset(mat,mat[,2]==subMat2)
	 									}
	 							}
						}
					endNodeL = rowSums(ll)
		 			mat = cbind(mat, endNodeL) # mat is now the full matrix with the root to node distance for each node
					startNodeL = matrix(1:l,nrow=l,ncol=1)
					startNodeL[] = 0
					for (j in 1:l)
						{
							r = mat[j,1]
							s = subset(mat,mat[,2]==r)
							for (k in 1:l)
								{
									if (nrow(s) > 0)
										{
											startNodeL[j,1] = s[,8]
	 									}
	 							}	
		 				}
					mat = cbind(mat,startNodeL)
					maxEndLIndice = which.max(mat[,"endNodeL"])
	 				maxEndL = mat[maxEndLIndice,"endNodeL"]
	 				endYear = matrix(mat[,"endNodeL"]-maxEndL)
	 				endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
					startYear = matrix(mat[,"startNodeL"]-maxEndL)
	 				startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
	 				mat = cbind(mat,startYear,endYear)
	 				
	 				# (6) Adding surface distance and tree ID:
	 				# distances = matrix(0, nrow=dim(mat)[1], ncol=1)
	 				# colnames(distances) = c("greatCircleDist_km")
	 				treeIDs = matrix(0, nrow=dim(mat)[1], ncol=1)
	 				colnames(treeIDs) = c("treeID")
	 				for (j in 1:dim(mat)[1])
	 					{
	 						# x1 = cbind(mat[j,"startLon"],mat[j,"startLat"])
	 						# x2 = cbind(mat[j,"endLon"],mat[j,"endLat"])
	 						# distances[j] = rdist.earth(x1, x2, miles=F, R=NULL)
	 						treeIDs[j] = treeID
	 					}
	 				# mat = cbind(mat,distances,treeIDs)
	 				mat = cbind(mat,treeIDs)
	 		
	 				file.remove(properTreeName)
					file = as.character(paste(localTreesDirectory, "/TreeExtractions_", i, ".csv", sep=""))
					write.csv(mat, file, row.names=F, quote=F)
				}	else	{
					file.remove(properTreeName)
					print(paste("A problem occured with tree extraction ",i,sep=""))
				}
		}
}
