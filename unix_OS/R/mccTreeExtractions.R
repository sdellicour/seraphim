mccTreeExtractions = function(mcc_tre, mostRecentSamplingDatum)
	{
		mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=7)
		colnames(mcc_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat")
		mcc_tab[,c("node1","node2")] = mcc_tre$edge
		mcc_tab[,c("length")] = mcc_tre$edge.length
		for (i in 1:length(mcc_tre$annotations))
			{
				annotations = mcc_tre$annotations[[i]]
				mcc_tab[i,c("endLon","endLat")] = cbind(annotations$location2, annotations$location1)
			}
		for (i in 1:length(mcc_tre$annotations))
			{
				index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
				if (length(index) > 0)
					{
						mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
					}	else		{
						annotations = mcc_tre$root.annotation
						mcc_tab[i,c("startLon","startLat")] = cbind(annotations$location2, annotations$location1)
					}
			}
		l = length(mcc_tab[,1]); ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
		for (j in 1:l)
			{
				subMat = mcc_tab[j,2]
				subMat = subset(mcc_tab,mcc_tab[,2]==subMat)
				ll[j,1] = subMat[,3]
				subMat = subMat[1,1]
				subMat1 = subset(mcc_tab,mcc_tab[,2]==subMat)
				for (k in 1:l)
					{
						if (nrow(subMat1) > 0)
							{
								ll[j,k+1] = subMat1[,3]
	 							subMat2 = subMat1[1,1]
	 							subMat1 = subset(mcc_tab,mcc_tab[,2]==subMat2)
	 						}
	 				}
			}
		endNodeL = rowSums(ll)
		mcc_tab = cbind(mcc_tab, endNodeL)
		startNodeL = matrix(1:l,nrow=l,ncol=1)
		startNodeL[] = 0
		for (j in 1:l)
			{
				r = mcc_tab[j,1]
				s = subset(mcc_tab,mcc_tab[,2]==r)
				for (k in 1:l)
					{
						if (nrow(s) > 0)
							{
								startNodeL[j,1] = s[,8]
	 						}
	 				}	
			}
		colnames(startNodeL) = "startNodeL"
		mcc_tab = cbind(mcc_tab,startNodeL)
		maxEndLIndice = which.max(mcc_tab[,"endNodeL"])
		maxEndL = mcc_tab[maxEndLIndice,"endNodeL"]
	 	endYear = matrix(mcc_tab[,"endNodeL"]-maxEndL)
	 	endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
		startYear = matrix(mcc_tab[,"startNodeL"]-maxEndL)
	 	startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
	 	colnames(startYear) = "startYear"; colnames(endYear) = "endYear"
	 	mcc_tab = cbind(mcc_tab,startYear,endYear)
		mcc_tab = mcc_tab[order(mcc_tab[,"startYear"],decreasing=F),]
		mcc_tab1 = mcc_tab[1,]; mcc_tab2 = mcc_tab[2:dim(mcc_tab)[1],]
		mcc_tab2 = mcc_tab2[order(mcc_tab2[,"endYear"],decreasing=F),]
		mcc_tab = rbind(mcc_tab1, mcc_tab2)		
		tipLabels = matrix(nrow=dim(mcc_tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
		for (i in 1:length(mcc_tre$tip.label))
			{
				tipLabels[which(mcc_tab[,"node2"]==i),"tipLabel"] = gsub("'","",mcc_tre$tip.label[i])
			}
		mcc_tab = cbind(mcc_tab, tipLabels)
		return(mcc_tab)
	}
