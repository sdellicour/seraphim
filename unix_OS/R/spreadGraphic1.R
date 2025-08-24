spreadGraphic1 = function(mcc_tre, mcc_tab, includeRoot=T) {

	polygons = list(); c = 0; p = list(); n = 0
	if (includeRoot == TRUE)
		{
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_1`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_1`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_1`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_2`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_2`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_2`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_3`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_3`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_3`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_4`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_4`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_4`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_5`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_5`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_5`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_6`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_6`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_6`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_7`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_7`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_7`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_8`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_8`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_8`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_9`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_9`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_9`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			if (!is.null(mcc_tre$root.annotation$`location2_80%HPD_10`))
				{
					xS = unlist(mcc_tre$root.annotation$`location2_80%HPD_10`)
					yS = unlist(mcc_tre$root.annotation$`location1_80%HPD_10`)
					n = n+1; p[[n]] = Polygon(cbind(xS,yS))
				}
			ps = Polygons(p,1); sps = SpatialPolygons(list(ps))
			spdf = SpatialPolygonsDataFrame(sps, data.frame(ID=1:length(sps)))
			endYear = mcc_tab[which(!mcc_tab[,"node1"]%in%mcc_tab[,"node2"])[1],"startYear"]
			names(spdf) = round(endYear,3); c = c+1; polygons[[c]] = spdf
		}
	for (i in 1:length(mcc_tre$annotations))
		{
			if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_1`))
				{
					p = list()
					xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_1`)
					yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_1`)
					p[[1]] = Polygon(cbind(xS,yS)); n = 1
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_2`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_2`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_2`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_3`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_3`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_3`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_4`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_4`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_4`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_5`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_5`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_5`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_6`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_6`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_6`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_7`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_7`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_7`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_8`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_8`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_8`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_9`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_9`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_9`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}
					if (!is.null(mcc_tre$annotations[[i]]$`location2_80%HPD_10`))
						{
							xS = unlist(mcc_tre$annotations[[i]]$`location2_80%HPD_10`)
							yS = unlist(mcc_tre$annotations[[i]]$`location1_80%HPD_10`)
							n = n+1; p[[n]] = Polygon(cbind(xS,yS))
						}		
					ps = Polygons(p,1); sps = SpatialPolygons(list(ps))
					spdf = SpatialPolygonsDataFrame(sps, data.frame(ID=1:length(sps)))
					endYear = mcc_tab[which(mcc_tab[,"node2"]==mcc_tre$edge[i,2]),"endYear"]
					names(spdf) = round(endYear,3); c = c+1; polygons[[c]] = spdf
				}
		}
	years = rep(NA, length(polygons))
	for (i in 1:length(polygons)) years[i] = as.numeric(names(polygons[[i]]))
	polygons = polygons[order(years)]
	return(polygons)
}
