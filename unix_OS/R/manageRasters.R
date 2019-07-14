# List all available and recommended rasters
listAvailableRasters = function() {
    
}

# Download rasters and crop to extent specified. Extent: xmin, xmax, ymin, ymax
downloadRasters = function(ids=c(), extent=c()) {
    
}

# Download the categorical IGBP land cover and crop to extent specified. Extent: xmin, xmax, ymin, ymax
downloadLandCoverRaster=function(ids=c(), extent=c()) {
    
}

# Decrease rasters resolution by a factor 'R'
decreaseResolution = function(envVariables=list(), R=5) {
	envVariables1 = envVariables
	envVariables2 = list()
	for (h in 1:length(envVariables1))
		{ 
			rast1 = henvVariables1[[]]
			rast2 = rast1
			nrow(rast2) = nrow(rast1)/R
			ncol(rast2) = ncol(rast1)/R
			rast2[] = 0
			mat1 = as.matrix(rast1)
			mat2 = as.matrix(rast2)
			for (i in 1:dim(rast2)[1])
				{
					for (j in 1:dim(rast2)[2])
						{
							averagedValue = 0
							counter = 0
							nonNA = FALSE
							nberNA = 0
							for (k in 1:R)
								{
									for (l in 1:R)
										{
											x = ((i-1)*R) + k
											y = ((j-1)*R) + l
											if ((x<(dim(mat1)[1]+1)) & (y<(dim(mat1)[2]+1)))
												{
													counter = counter + 1	
													if (!is.na(mat1[x,y]))
														{
															averagedValue = averagedValue + mat1[x,y]
															# nonNA = TRUE
														}	else	{
															nberNA = nberNA + 1
														}	
												}
										}
								}
							if (nberNA > (counter/2))
								{
									mat2[i,j] = NA
								}	else	{
									mat2[i,j] = averagedValue/counter
								}
						}	
				}
			rast2[,] = mat2
			names(rast2) = names(rast1)
			envVariables2[[h]] = rast2
		}
	return(envVariables2)
}

# Scale rasters and return list of 'n' rasters for each scaling k parameter
scaleRasters = function(envVariables=list(), k=c(10,100,1000)) {
  	envVariables1 = envVariables
	envVariables2 = list(); c = 0
	for (h in 1:length(envVariables1))
		{
			for (i in 1:length(k))
				{
					r = henvVariables1[[h]]
					M = max(r[], na.rm=T); r[] = (r[]*(k[i]/M))+1
					names(r) = paste0(names(envVariables1[h]),"_k",k[i])
					c = c+1; envVariables2[[c]] = r
				}
		}
}

# Start from the IBGP land cover raster (MCD12Q1, LC1) and generate a distinct continuous raster for each selected land cover variable
prepareLandCoverRasters = function(landCoverRaster, R=5, variables=("forests","croplands","urban_areas")) {
	envVariables = list(); c = 0
	# Categorical IBGP values:
		# 0 - water
		# 1 - evergreen Needleleaf forest
		# 2 - evergreen Broadleaf forest
		# 3 - deciduous Needleleaf forest
		# 4 - deciduous Broadleaf forest
		# 5 - mixed forest
		# 6 - closed shrublands
		# 7 - open shrublands
		# 8 - woody savannas
		# 9 - savannas
		# 10 - grasslands
		# 11 - permanent wetlands
		# 12 - croplands
		# 13 - urban and built-up
		# 14 - cropland/Natural vegetation mosaic
		# 15 - snow and ice
		# 16 - barren or sparsely vegetated
		# 254 - unclassified
		# 255 - fill Value
	rast1 = landCoverRaster
	coverNames = list()
	coverNames[[1]] = "water"
	coverNames[[2]] = "forests"
	coverNames[[3]] = "shrublands"
	coverNames[[4]] = "savannas"
	coverNames[[5]] = "grasslands"
	coverNames[[6]] = "wetlands"
	coverNames[[7]] = "croplands"
	coverNames[[8]] = "urban_areas"
	coverNames[[9]] = "snow_ice"
	coverNames[[10]] = "barren_vegetation"
	rast1[rast1[]==0] = 10 # water
	rast1[rast1[]==1] = 20 # forests
	rast1[rast1[]==2] = 20 # forests
	rast1[rast1[]==3] = 20 # forests
	rast1[rast1[]==4] = 20 # forests
	rast1[rast1[]==5] = 20 # forests
	rast1[rast1[]==6] = 30 # shrublands
	rast1[rast1[]==7] = 30 # shrublands
	rast1[rast1[]==8] = 40 # savannas
	rast1[rast1[]==9] = 40 # savannas
	rast1[rast1[]==10] = 50 # grasslands
	rast1[rast1[]==11] = 60 # wetlands
	rast1[rast1[]==12] = 70 # croplands
	rast1[rast1[]==13] = 80 # urban areas
	rast1[rast1[]==14] = 70 # croplands
	rast1[rast1[]==15] = 90 # snow, ice
	rast1[rast1[]==16] = 100 # barren, sparse vegetation
	rast1[] = rast1[]/10
	indices = which(coverNames%in%variables)
	for (g in indices)
		{	
			rast2 = rast1
			nrow(rast2) = nrow(rast1)/R
			ncol(rast2) = ncol(rast1)/R
			rast2[] = 0
			if (g != 1)
				{
					mat1 = as.matrix(rast1)
				}	else	{
					temp1 = rast1
					temp1[is.na(temp1[])] = 1
					mat1 = as.matrix(temp1)
				}
			mat2 = as.matrix(rast2)
			for (i in 1:dim(rast2)[1])
				{
					for (j in 1:dim(rast2)[2])
						{
							averagedValue = 0
							counter = 0
							maxValue = 0
							nonNA = FALSE
							nberNA = 0
							for (k in 1:R)
								{
									for (l in 1:R)
										{
											x = ((i-1)*R) + k
											y = ((j-1)*R) + l
											if ((x<(dim(mat1)[1]+1)) & (y<(dim(mat1)[2]+1)))
												{
													counter = counter + 1	
													if (!is.na(mat1[x,y]))
														{
															if (mat1[x,y] == g)
																{
																	averagedValue = averagedValue + 1
																}
														}	else	{
															nberNA = nberNA + 1
														}			
												}
										}
								}
							if (nberNA > (counter/2))
								{
									mat2[i,j] = NA
								}	else	{
									averagedValue = (averagedValue/counter)*(R*R)	
									mat2[i,j] = averagedValue
								}
						}	
				}
			rast2[,] = mat2
    		if (g == 1) # water  
    			{
					rast2[is.na(rast2)] = max(rast2[], na.rm=T)
      			}
			names(rast2) = paste("Land_cover",coverNames[[g]],suffix,sep="_")
			fileName = paste(names(rast2), ".asc", sep="")
			c = c+1; envVariables[[c]] = rast2
		}
}
