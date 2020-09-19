torusRandomisation <-
function(rast)	{

	nberOfRows = dim(rast)[1]; nberOfColumns = dim(rast)[2]
	jumpOfRows = as.integer(runif (1, min = 1, nberOfRows))
	jumpOfColumns = as.integer(runif (1, min = 1, nberOfColumns))
	reflectLines = FALSE; random = runif(1)
	if (random < 0.5)
		{
			reflectLines = TRUE
		}
	reflectColumns = FALSE; random = runif(1)
	if (random < 0.5)
		{	
			reflectColumns = TRUE
		}						
	mat = matrix(0, nrow=dim(rast)[1], ncol=dim(rast)[2])
	mat_transit = matrix(0, nrow=dim(rast)[1], ncol=dim(rast)[2])
	mat = raster::as.matrix(rast)
	mat_transit = mat
	random = runif(1); randomTot = FALSE
	randomColumns = FALSE; randomRows = FALSE
	counter = 0
	while (randomTot == FALSE)
		{	
			if ((random < 0.5) & (randomColumns == FALSE))
				{
					counter = counter + 1
					for (i in 1:nberOfRows)
						{
							NAs = 0
							for (j in 1:nberOfColumns)
								{
									if (!is.na(mat[i,j]))
										{
											newJ = j + jumpOfColumns + NAs
											if ((j + jumpOfColumns + NAs) > nberOfColumns)
												{
													newJ = j + jumpOfColumns + NAs - nberOfColumns
												}
											while (is.na(mat_transit[i,newJ]))
												{
													NAs = NAs + 1
													newJ = newJ + 1
													if (newJ > nberOfColumns)
														{
															newJ = newJ - nberOfColumns
														}
													if (newJ < 0)
														{
															newJ = newJ + nberOfColumns
														}
												}
											newJ = j + jumpOfColumns + NAs
											if ((j + jumpOfColumns + NAs) > nberOfColumns)
												{ 
													newJ = j + jumpOfColumns + NAs - nberOfColumns
												}
											if ((j + jumpOfColumns + NAs) < 0)
												{ 
													newJ = j + jumpOfColumns + NAs + nberOfColumns
												}
											mat[i,j] = mat_transit[i,newJ]
										}	else	{
											NAs = NAs - 1
										}
								}
							if (reflectLines == TRUE)
								{
									NAs = 0	
									line_transit = matrix(nrow=1, ncol=dim(rast)[2])
									for (j in 1:nberOfColumns)
										{
											if (!is.na(mat[i,j]))
												{
													J = nberOfColumns - (j-1) - NAs
													while(is.na(mat[i,J]))
														{
															NAs = NAs + 1
															J = nberOfColumns - (j-1) - NAs
														}
													line_transit[1,j] = mat[i,J]
												}	else	{
													NAs = NAs - 1
												}
										}
									mat[i,] = line_transit
								}		
						}
					mat_transit = mat
					randomColumns = TRUE
					random = 1
				}	else if (randomRows == FALSE)	{
					counter = counter + 1											
					for (j in 1:nberOfColumns)
						{
							NAs = 0
							for (i in 1:nberOfRows)
								{
									if (!is.na(mat[i,j]))
										{
											newI = i + jumpOfRows + NAs
											if ((i + jumpOfRows + NAs) > nberOfRows)
												{
													newI = i + jumpOfRows + NAs - nberOfRows
												}
											while (is.na(mat_transit[newI,j]))
												{
													NAs = NAs + 1
													newI = newI + 1
													if (newI > nberOfRows)
														{
															newI = newI - nberOfRows
														}
													if (newI < 0)
														{
															newI = newI + nberOfRows
														}
												}
											newI = i + jumpOfRows + NAs
											if ((i + jumpOfRows + NAs) > nberOfRows)
												{ 
													newI = i + jumpOfRows + NAs - nberOfRows
												}
											if ((i + jumpOfRows + NAs) < 0)
												{ 
													newI = i + jumpOfRows + NAs + nberOfRows
												}	
											mat[i,j] = mat_transit[newI,j]
										}	else	{
											NAs = NAs - 1
										}
								}
							if (reflectColumns == TRUE)
								{
									NAs = 0	
									column_transit = matrix(nrow=dim(rast)[1], ncol=1)
									for (i in 1:nberOfRows)
										{
											if (!is.na(mat[i,j]))
												{
													I = nberOfRows - (i-1) - NAs
													while(is.na(mat[I,j]))
														{
															NAs = NAs + 1
															I = nberOfRows - (i-1) - NAs
														}
													column_transit[i,1] = mat[I,j]
												}	else	{
													NAs = NAs - 1
												}
										}
									mat[,j] = column_transit
								}	
						}
					mat_transit = mat
					randomRows = TRUE
					random = 0
				}
			if (randomColumns == TRUE)
				{
					if (randomRows == TRUE)
						{
							randomTot = TRUE
						}
				}
		}
	rast[,] = mat
	return(rast)
}
