variogramModel <-
function(envVariables)\t{

\tvariogramModels = list()
\tfor (h in 1:length(envVariables))
\t\t{
\t\t\tpath = envVariables[[h]]@data@names
\t\t\tpath = unlist(strsplit(path, "/"))
\t\t\tpath = path[length(path)]
\t\t\tpath = unlist(strsplit(path, "\\\\."))
\t\t\tpath = gsub(" ", "_", path[1])
\t\t\tnames(envVariables[[h]]) = path
\t\t}
\tnscore = function(x)
\t\t{
\t\t\t# Function created by Ashton Shortridge (May/June 2008).\t
\t   \t\t# Takes a vector of values x and calculates their normal scores. Returns a list with the scores
\t   \t\t# and an ordered table of original values and scores, which is useful as a back-transform table.
\t  \t\tnscore = qqnorm(x, plot.it=F)$x  # normal score 
\t   \t\ttrn.table = data.frame(x=sort(x), nscore=sort(nscore))
\t   \t\treturn(list(nscore=nscore, trn.table=trn.table))
\t\t}
\tnscoreBack = function(scores, nscore)
\t\t{
\t\t\t# Function created by Ashton Shortridge (May/June 2008).\t\t
\t    \t# Given a vector of normal scores and a normal score object (from nscore), the function returns
\t   \t\t# a vector of back-transformed values.
\t\t\tmin.x = min(nscore$trn.table$x, na.rm=T)
\t \t\tmax.x = max(nscore$trn.table$x, na.rm=T)
\t   \t\tmin.sc = min(scores, na.rm=T)
\t  \t\tmax.sc = max(scores, na.rm=T)
\t\t\tx = c(min.x, nscore$trn.table$x, max.x)
\t    \tnsc = c(min.sc, nscore$trn.table$nscore, max.sc)
\t\t\tback.xf = approxfun(nsc,x) # Develop the back transform function
\t   \t\tval = back.xf(scores)
\t\t\treturn(val)
\t\t}
\tfor (i in 1:length(envVariables))
\t\t{
\t\t\t# 1. Sample the raster with 10,000 points:
\t\t\t\trast = envVariables[[i]]; # plot(rast, useRaster = F)
\t\t\t\t
\t\t\t# 2. Normal score transformation:
\t\t\t\trastScore = nscore(rast[]); rast2 = rast; rast2[] = rastScore$nscore
\t\t\t\t# rast2[] = nscoreBack(rast2[])
\t\t\t\trast = rast2 # plot(rast2)
\t\t\t\t
\t\t\t# 3. Compute a semivariogram of this sample:
\t\t\t\tif (ncell(rast) > 10^4)
\t\t\t\t\t{
\t\t\t\t\t\tnberOfPoints = 10^4
\t\t\t\t\t}\telse\t{
\t\t\t\t\t\tnberOfPoints = ncell(rast)
\t\t\t\t\t}\t
\t\t\t\tsamples = as.data.frame(sampleRandom(rast, nberOfPoints, sp=T))
\t\t\t\tcoordinates(samples) = ~x+y 
\t\t\t\t# also transforms into a "SpatialPointsDataFrame"
\t\t\t\t# plot(samples, add = T, pch=15, cex=0.2)
\t\t\t\tfmla = as.formula(paste(rast@data@names, "~1", sep=""))
\t\t\t\tobservedVariogram = variogram(fmla, data = samples)
\t\t\t\tdev.new(width=5, height =5); par("mar"=c(3.5,3.5,4.5,2)); par(mgp=c(1.6, 0.5, 0))
\t\t\t\ttext1 = paste("Observed variogram for ", rast@data@names, "*", sep="")
\t\t\t\ttext2 = paste("User-defined variogram for ", rast@data@names, "*", sep="")
\t\t\t\ttext3 = paste("Fitted (corrected) variogram for ", rast@data@names, "*", sep="")
\t\t\t\ttext4 = paste("* a normal transformation was                 ", sep="")
\t\t\t\ttext5 = paste(" previously performed on the raster grid   ", sep="")
\t\t\t\tplot(observedVariogram$dist, observedVariogram$gamma, xlab="distance", ylab="semivariance", cex.lab=0.7, cex.axis=0.7)
\t\t\t\tmtext(text1, cex=0.8, line=2.5)
\t\t\t\tmtext(text4, cex=0.7, line=-16, adj=1)
\t\t\t\tmtext(text5, cex=0.7, line=-16.7, adj=1)
\t\t\t\thappyWithModel = FALSE
\t\t\t\twhile(happyWithModel == FALSE)
\t\t\t\t\t{
\t\t\t\t\t\tcat("Setting parameters value for variogram model fitting (", rast@data@names, "):", sep="")
\t\t\t\t\t\tcat ("N.B: call vgm() without a model argument to get available models");
\t\t\t\t\t\tmodel = readline(prompt = "\tmodel (e.g. Exp, Sph, Gau, Mat) = ")
\t\t\t\t\t\tsill = readline(prompt = "\tsill value = ")
\t\t\t\t\t\tsill = as.numeric(sill)
\t\t\t\t\t\trange = readline(prompt = "\trange value = ")
\t\t\t\t\t\trange = as.numeric(range)
\t\t\t\t\t\tnugget = readline(prompt = "\tnugget value = ")
\t\t\t\t\t\tnugget = as.numeric(nugget)
\t\t\t\t\t\tuserVariogramModel = vgm(sill, model, range, nugget)
\t\t\t\t\t\tfittedVariogramModel = fit.variogram(observedVariogram, userVariogramModel)
\t\t\t\t\t\t\t# vgm: (sill, model, range, nugget)
\t\t\t\t\t\t\t# fit.variogram: fit sill and range
\t\t\t\t\t\t# plot.Variogram(observedVariogram, userVariogramModel)
\t\t\t\t\t\t# if (i != length(envVariables))
\t\t\t\t\t\t#\t{
\t\t\t\t\t\t#\t\tcat ("Press ENTER to continue"); line=readline()
\t\t\t\t\t\t#\t}
\t\t\t\t\t\tmtext(text2, cex=0.8, line=1.6, col='blue')
\t\t\t\t\t\tlines(variogramLine(userVariogramModel, max(observedVariogram$dist)), col='blue', lwd=2, lty=2)
\t\t\t\t\t\tSys.sleep(2)
\t\t\t\t\t\tmtext(text3, cex=0.8, line=0.7, col='red')
\t\t\t\t\t\tlines(variogramLine(fittedVariogramModel, max(observedVariogram$dist)), col='red', lwd=2)
\t\t\t\t\t\tyesOrNo = readline(prompt = "Are you satisfayed by the fitted variogram model in red (y/n)? ")
\t\t\t\t\t\tif (yesOrNo == "y")
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\thappyWithModel = TRUE
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tdev.off()
\t\t\t# 4. Create a mask where to store simulated data:
\t\t\t\t# mask = rast*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
\t\t\t# 5. Run the sequential simulation with the same range, sill and nugget as temperature:
\t\t\t\tvariogramModel = gstat(formula=fmla, dummy=T, beta=0, model=fittedVariogramModel, nmax=round(fittedVariogramModel$range[2]*10))
\t\t\t\tvariogramModels[[i]] = variogramModel
\t\t\t\tnScoreTransformation = FALSE; quantileTransformation = FALSE
\t\t\t\t# 1st OPTION: normal score transformation
\t\t\t\tif(nScoreTransformation == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\trastScore = nscore(rast[]); simRast = rast; simRast[] = rastScore$nscore
\t\t\t\t\t\tmask = simRast*0+1; mask.grid = as(mask, 'SpatialGridDataFrame')
\t\t\t\t\t\tsimulation = predict(variogramModels[[h]], newdata=mask.grid, nsim=1)
\t\t\t\t\t\tsimRast = raster(simulation)*mask
\t\t\t\t\t\tsimRast[] = nscoreBack(simRast[],rastScore)
\t\t\t\t\t}
\t\t\t\t# 2nd OPTION: quantile trasnformation
\t\t\t\tif(nScoreTransformation == TRUE)
\t\t\t\t\t{\t\t\t\t\t\t
\t\t\t\t\t\tsimulation = predict(variogramModel, newdata=mask.grid, nsim=1)
\t\t\t\t\t\tsimRast = raster(simulation)*mask #\tplot(simRast)
\t\t\t\t\t\trastData = as.data.frame(rast)
\t\t\t\t\t\tnames(rastData) = "value"
\t\t\t\t\t\trastData = rastData[with(rastData, order(value)),]
\t\t\t\t\t\tsimRastData = as.data.frame(simRast)
\t\t\t\t\t\tnames(simRastData) = "value"
\t\t\t\t\t\tindexes = 1:dim(simRastData)[1]
\t\t\t\t\t\tsimRastData$indexes = indexes
\t\t\t\t\t\tsimRastData = simRastData[with(simRastData, order(value)),]
\t\t\t\t\t\tsimRastData$value = rastData
\t\t\t\t\t\tsimRastData = simRastData[with(simRastData, order(indexes)),]
\t\t\t\t\t\tsimRast[] = simRastData[,"value"] #\tplot(simRast)
\t\t\t\t\t\trangeSimRast = max(simRast[], na.rm=T) - min(simRast[], na.rm=T)
\t\t\t\t\t\trangeRast = max(rast[], na.rm=T) - min(rast[], na.rm=T)
\t\t\t\t\t\tsimRast[] = simRast[] - min(simRast[], na.rm=T)
\t\t\t\t\t\tsimRast[] = (simRast[]/rangeSimRast)*rangeRast
\t\t\t\t\t\tsimRast[] = simRast[] + min(rast[], na.rm=T)
\t\t\t\t\t}
\t\t\t\t# Graphical comparison:
\t\t\t\t\t# plot(rast)
\t\t\t\t\t# plot(simRast)
\t\t\t\t\t# hist(rast)
\t\t\t\t\t# hist(simRast)
\t\t\t\t\t# Moran(rast)
\t\t\t\t\t# Moran(simRast)
\t\t\t\t\t# plot(MoranLocal(rast))
\t\t\t\t\t# plot(MoranLocal(simRast))
\t\t\t\t\t# Geary(rast)
\t\t\t\t\t# Geary(simRast)
\t\t\t\t\t# plot(GearyLocal(rast))
\t\t\t\t\t# plot(GearyLocal(simRast))
\t\t}
\treturn(variogramModels)
}
