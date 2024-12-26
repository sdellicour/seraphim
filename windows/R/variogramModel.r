variogramModel = function(envVariables)	{

	variogramModels = list()
	for (h in 1:length(envVariables))
		{
			path = envVariables[[h]]@data@names
			path = unlist(strsplit(path, "/"))
			path = path[length(path)]
			path = unlist(strsplit(path, "\\."))
			path = gsub(" ", "_", path[1])
			names(envVariables[[h]]) = path
		}
	nscore = function(x)
		{
			# Function created by Ashton Shortridge (May/June 2008).	
	   		# Takes a vector of values x and calculates their normal scores. Returns a list with the scores
	   		# and an ordered table of original values and scores, which is useful as a back-transform table.
	  		nscore = qqnorm(x, plot.it=F)$x  # normal score 
	   		trn.table = data.frame(x=sort(x), nscore=sort(nscore))
	   		return(list(nscore=nscore, trn.table=trn.table))
		}
	nscoreBack = function(scores, nscore)
		{
			# Function created by Ashton Shortridge (May/June 2008).		
	    	# Given a vector of normal scores and a normal score object (from nscore), the function returns
	   		# a vector of back-transformed values.
			min.x = min(nscore$trn.table$x, na.rm=T)
	 		max.x = max(nscore$trn.table$x, na.rm=T)
	   		min.sc = min(scores, na.rm=T)
	  		max.sc = max(scores, na.rm=T)
			x = c(min.x, nscore$trn.table$x, max.x)
	    	nsc = c(min.sc, nscore$trn.table$nscore, max.sc)
			back.xf = approxfun(nsc,x) # Develop the back transform function
	   		val = back.xf(scores)
			return(val)
		}
	for (i in 1:length(envVariables))
		{
			# 1. Sample the raster with 10,000 points:
				rast = envVariables[[i]]; # plot(rast, useRaster = F)
				
			# 2. Normal score transformation:
				rastScore = nscore(rast[]); rast2 = rast; rast2[] = rastScore$nscore
				# rast2[] = nscoreBack(rast2[])
				rast = rast2 # plot(rast2)
				
			# 3. Compute a semivariogram of this sample:
				if (ncell(rast) > 10^4)
					{
						nberOfPoints = 10^4
					}	else	{
						nberOfPoints = ncell(rast)
					}	
				samples = as.data.frame(sampleRandom(rast, nberOfPoints, sp=T))
				coordinates(samples) = ~x+y 
				# also transforms into a "SpatialPointsDataFrame"
				# plot(samples, add = T, pch=15, cex=0.2)
				fmla = as.formula(paste(rast@data@names,"~1",sep=""))
				observedVariogram = variogram(fmla, data = samples)
				dev.new(width=5, height =5); par("mar"=c(3.5,3.5,4.5,2)); par(mgp=c(1.6, 0.5, 0))
				text1 = paste("Observed variogram for ",rast@data@names,"*",sep="")
				text2 = paste("User-defined variogram for ",rast@data@names,"*",sep="")
				text3 = paste("Fitted (corrected) variogram for ",rast@data@names,"*",sep="")
				text4 = paste("* a normal transformation was                 ",sep="")
				text5 = paste(" previously performed on the raster grid   ",sep="")
				plot(observedVariogram$dist, observedVariogram$gamma, xlab="distance", ylab="semivariance", cex.lab=0.7, cex.axis=0.7)
				mtext(text1, cex=0.8, line=2.5)
				mtext(text4, cex=0.7, line=-16, adj=1)
				mtext(text5, cex=0.7, line=-16.7, adj=1)
				happyWithModel = FALSE
				while(happyWithModel == FALSE)
					{
						cat("Setting parameters value for variogram model fitting (", rast@data@names, "):",sep="")
						cat ("N.B: call vgm() without a model argument to get available models");
						model = readline(prompt = "	model (e.g. Exp, Sph, Gau, Mat) = ")
						sill = readline(prompt = "	sill value = ")
						sill = as.numeric(sill)
						range = readline(prompt = "	range value = ")
						range = as.numeric(range)
						nugget = readline(prompt = "	nugget value = ")
						nugget = as.numeric(nugget)
						userVariogramModel = vgm(sill, model, range, nugget)
						fittedVariogramModel = fit.variogram(observedVariogram, userVariogramModel)
							# vgm: (sill, model, range, nugget)
							# fit.variogram: fit sill and range
						# plot.Variogram(observedVariogram, userVariogramModel)
						# if (i != length(envVariables))
						#	{
						#		cat ("Press ENTER to continue"); line=readline()
						#	}
						mtext(text2, cex=0.8, line=1.6, col='blue')
						lines(variogramLine(userVariogramModel, max(observedVariogram$dist)), col='blue', lwd=2, lty=2)
						Sys.sleep(2)
						mtext(text3, cex=0.8, line=0.7, col='red')
						lines(variogramLine(fittedVariogramModel, max(observedVariogram$dist)), col='red', lwd=2)
						yesOrNo = readline(prompt = "Are you satisfayed by the fitted variogram model in red (y/n)? ")
						if (yesOrNo == "y")
							{
								happyWithModel = TRUE
							}
					}
				dev.off()
			# 4. Create a mask where to store simulated data:
				# mask = rast*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
			# 5. Run the sequential simulation with the same range, sill and nugget as temperature:
				variogramModel = gstat(formula=fmla, dummy=T, beta=0, model=fittedVariogramModel, nmax=round(fittedVariogramModel$range[2]*10))
				variogramModels[[i]] = variogramModel
				nScoreTransformation = FALSE; quantileTransformation = FALSE
				# 1st OPTION: normal score transformation
				if(nScoreTransformation == TRUE)
					{
						rastScore = nscore(rast[]); simRast = rast; simRast[] = rastScore$nscore
						mask = simRast*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
						simulation = predict(variogramModels[[h]], newdata=mask.grid, nsim=1)
						simRast = raster(simulation)*mask
						simRast[] = nscoreBack(simRast[],rastScore)
					}
				# 2nd OPTION: quantile trasnformation
				if(nScoreTransformation == TRUE)
					{						
						simulation = predict(variogramModel, newdata=mask.grid, nsim=1)
						simRast = raster(simulation)*mask #	plot(simRast)
						rastData = as.data.frame(rast)
						names(rastData) = "value"
						rastData = rastData[with(rastData, order(value)),]
						simRastData = as.data.frame(simRast)
						names(simRastData) = "value"
						indexes = 1:dim(simRastData)[1]
						simRastData$indexes = indexes
						simRastData = simRastData[with(simRastData, order(value)),]
						simRastData$value = rastData
						simRastData = simRastData[with(simRastData, order(indexes)),]
						simRast[] = simRastData[,"value"] #	plot(simRast)
						rangeSimRast = max(simRast[], na.rm=T) - min(simRast[], na.rm=T)
						rangeRast = max(rast[], na.rm=T) - min(rast[], na.rm=T)
						simRast[] = simRast[] - min(simRast[], na.rm=T)
						simRast[] = (simRast[]/rangeSimRast)*rangeRast
						simRast[] = simRast[] + min(rast[], na.rm=T)
					}
				# Graphical comparison:
					# plot(rast)
					# plot(simRast)
					# hist(rast)
					# hist(simRast)
					# Moran(rast)
					# Moran(simRast)
					# plot(MoranLocal(rast))
					# plot(MoranLocal(simRast))
					# Geary(rast)
					# Geary(simRast)
					# plot(GearyLocal(rast))
					# plot(GearyLocal(simRast))
		}
	return(variogramModels)
}    