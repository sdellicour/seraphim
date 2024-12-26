postTreeExtractions = function(post_tre, mostRecentSamplingDatum)
	{
		post_tab = matrix(nrow=dim(post_tre$edge)[1], ncol=7)
		colnames(post_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat")
		post_tab[,c("node1","node2")] = post_tre$edge
		post_tab[,c("length")] = post_tre$edge.length
		for (i in 1:length(post_tre$annotations))
			{
				annotations = post_tre$annotations[[i]]
				post_tab[i,c("endLon","endLat")] = cbind(annotations$location[[2]], annotations$location[[1]])
			}
		for (i in 1:length(post_tre$annotations))
			{
				index = which(post_tab[,"node2"] == post_tab[i,"node1"])
				if (length(index) > 0)
					{
						post_tab[i,c("startLon","startLat")] = post_tab[index,c("endLon","endLat")]
					}	else		{
						annotations = post_tre$root.annotation
						post_tab[i,c("startLon","startLat")] = cbind(annotations$location[[2]], annotations$location[[1]])
					}
			}
		l = length(post_tab[,1]); ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
		for (j in 1:l)
			{
				subMat = post_tab[j,2]
				subMat = subset(post_tab,post_tab[,2]==subMat)
				ll[j,1] = subMat[,3]
				subMat = subMat[1,1]
				subMat1 = subset(post_tab,post_tab[,2]==subMat)
				for (k in 1:l)
					{
						if (nrow(subMat1) > 0)
							{
								ll[j,k+1] = subMat1[,3]
	 							subMat2 = subMat1[1,1]
	 							subMat1 = subset(post_tab,post_tab[,2]==subMat2)
	 						}
	 				}
			}
		endNodeL = rowSums(ll)
		post_tab = cbind(post_tab, endNodeL)
		startNodeL = matrix(1:l,nrow=l,ncol=1)
		startNodeL[] = 0
		for (j in 1:l)
			{
				r = post_tab[j,1]
				s = subset(post_tab,post_tab[,2]==r)
				for (k in 1:l)
					{
						if (nrow(s) > 0)
							{
								startNodeL[j,1] = s[,8]
	 						}
	 				}	
			}
		colnames(startNodeL) = "startNodeL"
		post_tab = cbind(post_tab,startNodeL)
		maxEndLIndice = which.max(post_tab[,"endNodeL"])
		maxEndL = post_tab[maxEndLIndice,"endNodeL"]
	 	endYear = matrix(post_tab[,"endNodeL"]-maxEndL)
	 	endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
		startYear = matrix(post_tab[,"startNodeL"]-maxEndL)
	 	startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
	 	colnames(startYear) = "startYear"; colnames(endYear) = "endYear"
	 	post_tab = cbind(post_tab,startYear,endYear)
		post_tab = post_tab[order(post_tab[,"startYear"],decreasing=F),]
		post_tab1 = post_tab[1,]; post_tab2 = post_tab[2:dim(post_tab)[1],]
		post_tab2 = post_tab2[order(post_tab2[,"endYear"],decreasing=F),]
		post_tab = rbind(post_tab1, post_tab2)		
		tipLabels = matrix(nrow=dim(post_tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
		for (i in 1:length(post_tre$tip.label))
			{
				tipLabels[which(post_tab[,"node2"]==i),"tipLabel"] = gsub("'","",post_tre$tip.label[i])
			}
		post_tab = cbind(post_tab, tipLabels)
		return(post_tab)
	}
