# This R script described how to generate the three panels of Figure 1 included in the application note of "seraphim" 2.0.
# Figure 1 constitutes an example of visualisations that can be generated with the package. These Visualisations are based
# on a continuous phylogeographic analysis of the yellow fever virus (YFV) outbreak that started around 2015 in
# southeastern Brazil (Hill et al. 2022, biorXiv). Figure 1A displays a continuous phylogeographic reconstruction of the
# dispersal history of YFV outbreak lineages, i.e. the maximum clade credibility (MCC) tree and overall 80% highest
# posterior density (HPD) regions reflecting the uncertainty of the Bayesian phylogeographic inference summarised from
# 1,000 trees sampled from the post-burn-in posterior tree distribution. MCC tree nodes are coloured according to their
# time of occurrence and 80% HPD regions were computed for successive time layers and then superimposed using the same
# colour scale to reflect time. Figure 1B shows the evolution of the maximal wavefront distance from the epidemic origin:
# the solid curve represents the median value and the surrounding polygon the 95% HPD interval. Those estimates are also
# based on 1,000 trees sampled from the post-burn-in posterior tree distribution. Finally, Figure 1C reports the
# evaluation of the diffusion velocity of viral lineages through the estimation of the weighted diffusion coefficient
# (WDC): kernel density estimates of the diffusion coefficient (DC) parameters, with the posterior WDC estimates on the
# x-axis and the coefficient of variation of the diffusion coefficient among the branches of each sampled tree on the
# y-axis (the shades of decreasing darkness showing the 50%, 75%, and 95% HPD regions via kernel density estimation).


library(seraphim) # to load the "seraphim" package, which can be installed from its GitHub repository
library(diagram) # additional package used for mapping phylogenetic branches as curves with the function "curvedarrow"


# Step 1: loading shapefiles withe Brazilian international and municipality borders on the study area (source: GADM)

study_area_country_border = shapefile("Shapefiles_Brazil/Study_area_country_border.shp")
study_area_municipalities = shapefile("Shapefiles_Brazil/Study_area_municipalities.shp")
extent_study_area = extent(study_area_municipalities) # to retrieve the extent of the study area


# Step 2: extracting the spatio-temporal information embedded in the maximum clade credibility (MCC) tree

mostRecentSamplingDatum = 2019.288
tree = readAnnotatedNexus("YFV_MCCtree.tree")
mcc = mccTreeExtractions(tree, mostRecentSamplingDatum)
write.csv(mcc, "YFV_MCCtree.csv", row.names=F, quote=F)


# Step 3: estimating the HPD regions for successive time slices to reflect the uncertainty of the Bayesian inference

localTreesDirectory = "Tree_extractions" # directory where the tree extraction tables are stored
nberOfExtractionFiles = 1000 # number of tree extraction tables gathered in that directory
prob = 0.80 # probability corresponding to the HPD (highest posterior density) regions.
precision = 0.5 # time interval that will be used to define the successive time slices
startDate = 2014.305 # corresponds to lower limit of the 95% HPD interval estimated for the root age
hpd_regions = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDate, precision))


# Step 4: defining the different colour scales and colour codes to use for the tree nodes and HPD regions

colour_scale = paste0(colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:111],"CC") # overall colour scale retrived
		# from the "RColorBrewer" package; "CC" allowing to define a level of transparency set to 20%
mcc = read.csv("YFV_MCCtree.csv", head=T) # to reload the tree extraction file dataframe generated in step #2
minYear = startDate; maxYear = max(mcc[,"endYear"]) # to retrieve the minimum and maximal decimal years
nodes_colours = colour_scale[(((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1] # assigns a colour to each node
root_node_index = (((mcc[which(!mcc[,"node1"]%in%mcc[,"node2"])[1],"startYear"]-minYear)/(maxYear-minYear))*100)+1
root_node_colour = colour_scale[root_node_index] # assigns a colour code to the root node of the tree
hpd_region_colours = rep(NA, length(hpd_regions)) # vector that will contain the colour code of each HPD region
for (i in 1:length(hpd_regions))
	{
		date = as.numeric(names(hpd_regions[[i]])) # date of HPD region i
		hpd_region_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		hpd_region_colours[i] = gsub("CC","4D",colour_scale[hpd_region_index])
				# "4D" allowing to define a level of transparency set to 70%
	}


# Step 5: plotting the dispersal history of YFV lineages (and uncertainty) up to three successive points in time

pdf(paste0("Figure_1A_NEW.pdf"), width=9, height=3.6) # dev.new(width=9, height=3.6)
par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(0.5,1,0.75,1), mgp=c(0,0.4,0), lwd=0.4,
	bty="o", col="gray30", col.axis="gray30", fg="gray30")
cutOffs = c(2016.62, 2017.62, mostRecentSamplingDatum) # maximum dates considered to generate a snapshot by period
for (h in 1:length(cutOffs))
	{
		plot(study_area_country_border, border=NA, col="gray93"); CEX = 1.0
		plot(study_area_municipalities, border="white", col=NA, lwd=0.5, add=T)
		plot(study_area_country_border, border="gray50", col=NA, lwd=0.5, add=T)
		for (i in 1:length(hpd_regions)) # (1) to superimpose the different HPD region polygons on the study area
			{
				hpd_region_date = as.numeric(names(hpd_regions[[i]]))
				if (hpd_region_date <= cutOffs[h]) # to only plot the hpd region polygons up the maximum date
					{
						pol = hpd_regions[[i]]; crs(pol) = crs(study_area_country_border)
						plot(pol, axes=F, col=hpd_region_colours[i], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1]) # (2) to plot all the phylogeny branches as curved lines
			{
				if (mcc[i,"endYear"] <= cutOffs[h])  # to only plot the branches occurring up the maximum date
					{
						diagram::curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]),
											 cbind(mcc[i,"endLon"],mcc[i,"endLat"]),
											 arr.length=0, arr.width=0, lwd=0.3, lty=1,
											 lcol="gray30", arr.col=NA, arr.pos=F,
											 curve=0.1, dr=NA, endhead=F)
					}
			}
		for (i in dim(mcc)[1]:1) # (3) to plot all the tree nodes coloured according to time
			{
				if (i == 1) # to plot the root nodes
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col="gray70", cex=CEX)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=root_node_colour, cex=CEX)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.3, cex=CEX)
					}
				if (mcc[i,"endYear"] <= cutOffs[h]) # to plot all the other nodes occurring up the maximum date
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col="gray70", cex=CEX)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=nodes_colours[i], cex=CEX)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.3, cex=CEX)
					}
			}
		rect(extent_study_area@xmin, extent_study_area@ymin, extent_study_area@xmax, extent_study_area@ymax,
			 lwd=0.5, border="gray30", lty=1) # to add rectangular contour to the generated maps
		if (cutOffs[h] == mostRecentSamplingDatum) # to add a colour scale legend next to the final map
			{
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDate; rast[2] = max(mcc[,"endYear"])
				par(lwd=0.1, bty="o", col="gray30", col.axis="gray30", fg=NA)
				plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.420,0.900,0.130,0.145),
					 legend.args=list(text="", cex=0.7, line=0.4, col=NA), horizontal=T,
			 		 axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.4, tck=-0.8, col.axis=NA, col.tick=NA, line=0, mgp=c(0,0.12,0), at=seq(2014,2019,1)))
				par(lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
				plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.420,0.900,0.130,0.145),
					 legend.args=list(text="", cex=0.7, line=0.4, col="gray30"), horizontal=T,
			 		 axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.4, tck=-0.8, col.axis="gray30", col.tick="gray30", line=0, mgp=c(0,0.12,0), at=seq(2014,2019,1)))
			}
	}
dev.off()


# Step 6: estimating the dispersal statistics of YFV lineages with the "spreadStatistics" function

timeSlices = 100 # number of time slices that will be used to generate the maximal wavefront distance evolution plots
onlyTipBranches = FALSE # boolean variable defining if the statistic estimations have to be performed only while 						# considering the tip branches of the trees
showingPlots = FALSE # boolean variable specifying if the different plots have to be displayed or not.
outputName = paste0("Dispersal_stats/YFV") # name (prefix) to give to the different output files
nberOfCores = 1 # number of available cores to parallelise the computations (only work on Unix operating systems)
slidingWindow = 1/12 # sliding window, in units of time, that will be considered (here set to ~1 month)
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)


# Step 7: plotting the evolution through time of the maximal wavefront distance from the epidemic origin

pdf(paste0("Figure_1B_NEW.pdf"), width=6.5, height=2.5) # dev.new(width=6.5, height=2.5)
par(oma=c(0.5,1,0.75,1), mar=c(2,2,0.5,0.5), mgp=c(0,0.4,0), lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
swf_median = read.table(paste0("Dispersal_stats/YFV_median_spatial_wavefront_distance.txt"), header=T)
swf_95pHPD = read.table(paste0("Dispersal_stats/YFV_95%HPD_spatial_wavefront_distance.txt"), header=T)
swf = cbind(swf_median, swf_95pHPD[,2:3]); colnames(swf) = c("time","median","95pHDP_lower","95pHDP_upper")
swf = swf[which((swf[,"time"]>=minYear)&(swf[,"time"]<maxYear)),]
xMin = min(swf[,"time"]); xMax = max(swf[,"time"]); timeSlice = swf[2,"time"]-swf[1,"time"]; xMin = 2014.47; xMax = 2019.1
yMin = min(swf[,"95pHDP_lower"], na.rm=T); yMax = max(swf[,"95pHDP_upper"], na.rm=T); yMin = 0; yMax = 1480
colours = gsub("CC","CC",colour_scale[(((swf[,c("time")]-minYear)/(maxYear-minYear))*(length(colour_scale)-1))+1])
plot(swf[,"time"], swf[,"median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
xx_l = c(swf[,c("time")],rev(swf[,c("time")])); yy_l = c(swf[,"95pHDP_lower"],rev(swf[,"95pHDP_upper"]))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
for (i in 1:length(swf[,"time"]))
	{
		x1 = swf[i,"time"]-(timeSlice/2); x2 = swf[i,"time"]+(timeSlice/2)
		y1 = swf[i,"95pHDP_lower"]-250; y2 = swf[i,"95pHDP_upper"]+250
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="gray50", border=NA)
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colours[i], border=NA)
	}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
lines(swf[,"time"], swf[,"median"], lwd=1.0, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
for (i in c(2015:2019)) abline(v=i, lty=2, col="gray30", lwd=0.3)
axis(side=1, cex.axis=0.7, lwd=0.4, lwd.tick=0.4, tck=-0.030, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(0,0.10,0), at=seq(2014,2020))
axis(side=2, cex.axis=0.7, lwd=0.4, lwd.tick=0.4, tck=-0.030, col="gray30", col.axis="gray30", col.tick="gray30", mgp=c(1,0.28,0), at=seq(0,1600,400))
mtext(" Distance from epidemic origin (km)", side=2, col="gray30", cex=0.7, line=1.2, las=3)
dev.off()


# Step 8: evaluating of the diffusion velocity of viral lineages with the weighted diffusion coefficient

pdf(paste0("Figure_1C_NEW.pdf"), width=3.0, height=2.9) # dev.new(width=3.0, height=2.9)
par(oma=c(0,0,0,0), mar=c(2.5,2.4,0.5,0.5), lwd=0.4, bty="o", col="gray30", col.axis="gray30", fg="gray30")
three_colours = c("#FFFFFF",paste0(colorRampPalette(brewer.pal(7,"YlOrBr"))(7)[c(4:6)]))
three_colours = c("#FFFFFF",paste0(colorRampPalette(brewer.pal(7,"Reds"))(7)[c(2,4)]),"#DE4326")
three_colours = c("#FFFFFF",paste0(colorRampPalette(brewer.pal(7,"Reds"))(7)[c(4:6)]),"#DE4326")
statistics = read.table(paste0("Dispersal_stats/YFV_estimated_dispersal_statistics.txt"), header=T)
mean_WDCs = statistics[,5]/365.25 # to get mean WDC values in units of km2 per day
H = Hpi(cbind(mean_WDCs,statistics[,6])); kde = kde(cbind(mean_WDCs,statistics[,6]),H=H)
plot(kde, display="filled.contour2", cont=c(50,75,95), col=three_colours, axes=F, ann=F, xlim=c(22,44), ylim=c(3.5,23.5))
axis(side=1, lwd.tick=0.4, cex.axis=0.7, lwd=0, tck=-0.025, mgp=c(1,0.2,0), col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.4, cex.axis=0.7, lwd=0, tck=-0.025, mgp=c(1,0.3,0), col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab=expression("Weighted diffusion coefficient (km"^2*"/day)"), cex.lab=0.7, mgp=c(1.1,0,0), col.lab="gray30")
title(ylab="DC variation among branches", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
legend(cbind(37.2,24.0), c("95% HPD","75% HPD","50% HPD"), text.col="gray30", pch=16, pt.cex=1.4, col=three_colours[2:4], box.lty=0, cex=0.65, y.intersp=1.2)
legend(cbind(37.2,24.0), c("95% HPD","75% HPD","50% HPD"), text.col=NA, pch=1, pt.lwd=0.4, pt.cex=1.4, col="gray30", box.lty=0, cex=0.65, y.intersp=1.2)
dev.off()

