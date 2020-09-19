mdsTransformation <-
function(input, envVariables=list(), pathModel=2, resistances=c(T), avgResistances=c(T), fourCells=F, outputName="", OS="Unix")	{

	input0 = input; showingPlots = FALSE
	nberOfCores = 1; registerDoMC(cores=nberOfCores)
	textFile = FALSE; fastaAlignment = FALSE; treeFile = FALSE
	if (attr(input,"class") == "data.frame") textFile = TRUE
	if (attr(input,"class") == "DNAbin") fastaAlignment = TRUE
	if (attr(input,"class") == "phylo") treeFile = TRUE
	if (textFile == TRUE)
		{
			sequencesID = seq(1,dim(input)[1],1); years = rep(0,dim(input)[1]); locations1 = input[,2]; locations2 = input[,1]
			coordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
		}
	if (fastaAlignment == TRUE)
		{
			fastaChar = as.character(input)
			labels = rownames(fastaChar)
			names = matrix(nrow=length(labels), ncol=4)
			sequences = matrix(nrow=dim(fastaChar)[1], ncol=1)
			for (k in 1:dim(fastaChar)[1])
				{
					names[k,1:4] = unlist(strsplit(labels[k], "_"))
					sequence = ""
					for (l in 1:dim(fastaChar)[2]) sequence = paste(sequence,fastaChar[k,l],sep="")
					sequence = gsub("a","A",sequence); sequence = gsub("c","C",sequence)
					sequence = gsub("g","G",sequence); sequence = gsub("t","T",sequence)
					sequences[k,1] = sequence
				}
			sequencesID = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
			coordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
		}
	if (treeFile == TRUE)
		{
			tree = input
			labels = tree$tip.label
			labels = gsub("'","",labels)
			names = matrix(nrow=length(labels), ncol=4)
			for (k in 1:length(labels))
				{
					names[k,1:4] = unlist(strsplit(labels[k], "_"))
				}
			sequencesID = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
			coordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
		}
	nullRaster = envVariables[[1]]; nullRaster[!is.na(nullRaster[])] = 1
	names(nullRaster) = "nullRaster"
	buffer1 = list(nullRaster)
	buffer2 = list(T); buffer3 = list(T)
	for (i in 1:length(envVariables))
		{
			buffer1[[i+1]] = envVariables[[i]]
			buffer2[[i+1]] = resistances[[i]]
			buffer3[[i+1]] = avgResistances[[i]]
		}
	envVariables = list(); envVariables = buffer1
	resistances = list(); resistances = buffer2
	avgResistances = list(); avgResistances = buffer3
	rasterNames = c()
	for (i in 1:length(envVariables)) rasterNames = c(rasterNames, names(envVariables[[i]]))
	if (pathModel == 3)
		{
			if("CS_envVariables"%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), "CS_envVariables"))
			for (i in 1:length(envVariables))
				{
					# name = paste("CS_envVariables/", gsub("_0.","_0,",names(envVariables[[i]])), "_cs.asc", sep="")
					name = paste("CS_envVariables/", names(envVariables[[i]]), "_cs.asc", sep="")
					writeRaster(envVariables[[i]], name, overwrite=T)
				}
		}
	dists = list(); euclDis = matrix(0, nrow=dim(coordinates)[1], ncol=dim(coordinates)[1])
	buffer1 = coordinates[,1:2]; coordinates0 = coordinates
	coordinates = cbind(as.numeric(coordinates0[,3]),as.numeric(coordinates0[,4]))
	minEuclDis = sqrt(((coordinates[2,1]-coordinates[1,1])^2)+((coordinates[2,2]-coordinates[1,2])^2))
	maxEuclDis = sqrt(((coordinates[2,1]-coordinates[1,1])^2)+((coordinates[2,2]-coordinates[1,2])^2))
	for (j in 2:dim(coordinates)[1])
		{
			for (k in 1:j)
				{
					d = sqrt(((coordinates[j,1]-coordinates[k,1])^2)+((coordinates[j,2]-coordinates[k,2])^2))
					euclDis[k,j] = d; euclDis[j,k] = d
					if (minEuclDis > d) minEuclDis = d
					if (maxEuclDis < d) maxEuclDis = d
				}
		}
	hullRaster = nullRaster
	hull = chull(coordinates)
	hull = c(hull,hull[1])
	p = Polygon(coordinates[hull,])
	ps = Polygons(list(p),1)
	sps = SpatialPolygons(list(ps))
	pointsRaster = rasterize(coordinates, crop(hullRaster, sps, snap="out"))
	pointsRaster[!is.na(pointsRaster[])] = 0
	hullRaster = crop(hullRaster, sps, snap="out")
	bufferRaster = hullRaster
	hullRaster = mask(hullRaster, sps, snap="out")
	hullRaster[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
	buffer = list()
	buffer = foreach(i = 1:length(envVariables)) %dopar% {
	# for (i in 1:length(envVariables)) {
			sequences_to_keep = labels
			dist = matrix(0, nrow=length(coordinates[,1]), ncol=length(coordinates[,1]))
			if (pathModel == 1)
				{
					method = "SL"
					linesList = list()
					c = 0
					for (j in 2:length(coordinates[,1]))
						{
							for (k in 1:(j-1))
								{
									c = c+1
									points = rbind(coordinates[k,], coordinates[j,])
									linesList[[c]] = Lines(list(Line(points)),c)
								}
						}
					lines = SpatialLines(linesList)
					extractions = extract(envVariables[[i]], lines)
					c = 0
					for (j in 2:length(coordinates[,1]))
						{
							for (k in 1:(j-1))
								{
									c = c+1
									dist[k,j] = sum(extractions[[c]], na.rm=T)
									dist[j,k] = dist[k,j]
								}
						}
				}
			if (pathModel == 2)
				{
					method = "LC"
					if (resistances[[i]] == T)
						{
							trEnvVariable = transition(envVariables[[i]], function(x) 1/mean(x), directions=4)
						}	else	{
							trEnvVariable = transition(envVariables[[i]], mean, directions=4)
						}
					trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
					dist = costDistance(trEnvVariableCorr, as.matrix(coordinates), as.matrix(coordinates))
				}
			if (pathModel == 3)
				{
					method = "RW"				
					dist = circuitScape1(envVariables[[i]], paste("CS_envVariables/",names(envVariables[[i]]),"_cs",sep=""), 
										resistances[[i]], avgResistances[[i]], fourCells, coordinates, coordinates, OS, outputName, i)
					dist = as.matrix(dist)
				}
			dist[!is.finite(dist[])] = NA
			# To remove lines/columns of NA:
			lines_to_remove = c()
			for (j in 1:dim(dist)[1])
				{
					allNA = TRUE
					nas = 0
					for (k in 1:dim(dist)[2])
						{
							if (!is.na(dist[j,k])) 
								{
									if (dist[j,k]!=0) allNA = FALSE
								}	else	{
									nas = nas+1
								}
						}
					if (nas >= 3) allNA = TRUE
					if (allNA == TRUE)
						{
							lines_to_remove = c(lines_to_remove, j)
						}
				}
			if (length(lines_to_remove) > 0)
				{
					cat(paste("Removed sequences:","\\n",sep=""))
					for (j in 1:length(lines_to_remove))
						{
							cat(paste("   ", buffer1[j,1],"\\n",sep=""))	
						}
					dist = dist[-lines_to_remove,]; dist = dist[,-lines_to_remove]
					buffer2 = buffer1[-lines_to_remove,]
					sequences_to_keep = sequences_to_keep[-lines_to_remove]
				}	else	{
					buffer2 = buffer1
				}
			# To remove lines/columns with a few or a single NA value(s):
			lines_to_remove = c(); nberOfNA = c()
			for (j in 1:dim(dist)[1])
				{
					atLeastOneNA = FALSE; n = 0
					colIDs = c()
					for (k in 1:dim(dist)[2])
						{
							if (is.na(dist[j,k]))
								{
									atLeastOneNA = TRUE
									n = n+1
									colIDs = c(colIDs, k)
								}
						}
					# if (length(colIDs) > 0) print(colIDs)
					if (atLeastOneNA == TRUE)
						{
							lines_to_remove = c(lines_to_remove, j)
							nberOfNA = c(nberOfNA, n)
						}
				}
			if (length(lines_to_remove) > 0)
				{
					lines_to_remove = lines_to_remove[order(nberOfNA,decreasing=T)]
					stillNA = TRUE; J = 0
					while (stillNA == TRUE)
						{
							J = J+1
							test = dist[-lines_to_remove[1:J],-lines_to_remove[1:J]]
							atLeastOneNA = FALSE; n = 0
							for (j in 1:dim(test)[1])
								{
									for (k in 1:dim(test)[2])
										{
											if (is.na(test[j,k]))
												{
													atLeastOneNA = TRUE
													n = n+1
												}
										}
								}
							print(n)	
							if (atLeastOneNA == FALSE)
								{
									dist = dist[-lines_to_remove[1:J],]; dist = dist[,-lines_to_remove[1:J]]
									buffer3 = buffer2[-lines_to_remove[1:J],]
									sequences_to_keep = sequences_to_keep[-lines_to_remove]
									cat(paste("Removed sequences:","\\n",sep=""))
									for (j in 1:length(lines_to_remove[1:J]))
										{
											cat(paste("   ", buffer1[lines_to_remove[j],1],"\\n",sep=""))	
										}
									stillNA = FALSE
								}
						}
				}	else	{
					buffer3 = buffer2
				}
			minDist = min(dist, na.rm=T); maxDist = max(dist, na.rm=T)
			dist_mod = (((dist-minDist)/(maxDist-minDist))*(maxEuclDis-minEuclDis))+minEuclDis
			coordsMDS = cmdscale(dist_mod, eig=T, k=2); coordsMDS_mod = coordsMDS$points
			coordsMDS_mod[,1] = coordsMDS_mod[,1]-median(coordsMDS_mod[,1])
			coordsMDS_mod[,2] = coordsMDS_mod[,2]-median(coordsMDS_mod[,2])
			coordsMDS_mod1 = coordsMDS_mod; coordsMDS_mod2 = coordsMDS_mod
			coordsMDS_mod1[,1] = coordsMDS_mod1[,1]+median(coordinates[,1])
			coordsMDS_mod1[,2] = coordsMDS_mod1[,2]+median(coordinates[,2])
			coordsMDS_mod2 = coordsMDS_mod2[,2:1]
			coordsMDS_mod2[,1] = coordsMDS_mod2[,1]+median(coordinates[,1])
			coordsMDS_mod2[,2] = coordsMDS_mod2[,2]+median(coordinates[,2])
			coordsMDS_mod1a = coordsMDS_mod1; coordsMDS_mod1b = coordsMDS_mod1; coordsMDS_mod1c = coordsMDS_mod1; coordsMDS_mod1d = coordsMDS_mod1
			coordsMDS_mod2a = coordsMDS_mod2; coordsMDS_mod2b = coordsMDS_mod2; coordsMDS_mod2c = coordsMDS_mod2; coordsMDS_mod2d = coordsMDS_mod2
			coordsMDS_mod1b[,1] = coordsMDS_mod1b[,1]+(2*(mean(coordsMDS_mod1b[,1])-coordsMDS_mod1b[,1]))
			coordsMDS_mod1c[,2] = coordsMDS_mod1c[,2]+(2*(mean(coordsMDS_mod1c[,2])-coordsMDS_mod1c[,2]))
			coordsMDS_mod1d[,1] = coordsMDS_mod1d[,1]+(2*(mean(coordsMDS_mod1d[,1])-coordsMDS_mod1d[,1]))
			coordsMDS_mod1d[,2] = coordsMDS_mod1d[,2]+(2*(mean(coordsMDS_mod1d[,2])-coordsMDS_mod1d[,2]))
			coordsMDS_mod2b[,1] = coordsMDS_mod2b[,1]+(2*(mean(coordsMDS_mod2b[,1])-coordsMDS_mod2b[,1]))
			coordsMDS_mod2c[,2] = coordsMDS_mod2c[,2]+(2*(mean(coordsMDS_mod2c[,2])-coordsMDS_mod2c[,2]))
			coordsMDS_mod2d[,1] = coordsMDS_mod2d[,1]+(2*(mean(coordsMDS_mod2d[,1])-coordsMDS_mod2d[,1]))
			coordsMDS_mod2d[,2] = coordsMDS_mod2d[,2]+(2*(mean(coordsMDS_mod2d[,2])-coordsMDS_mod2d[,2]))
			coordsMDS_modList = list(coordsMDS_mod1a,coordsMDS_mod1b,coordsMDS_mod1c,coordsMDS_mod1d,coordsMDS_mod2a,coordsMDS_mod2b,coordsMDS_mod2c,coordsMDS_mod2d)
			numberOfNApoints = c()
			for (j in 1:length(coordsMDS_modList))
				{
					n = 0
					for (k in 1:dim(coordsMDS_modList[[j]])[1])
						{
							if(is.na(raster::extract(hullRaster,cbind(coordsMDS_modList[[j]][k,1],coordsMDS_modList[[j]][k,2])))) n = n+1
						}
					numberOfNApoints = c(numberOfNApoints, n)
				}
			if (showingPlots == TRUE)
				{	
					# plotRaser(nullRaster, col="grey90")
					# lines(coordinates[hull,], lwd=0.5); points(coordinates, pch=19, cex=0.7)
					# points(coordsMDS_mod1a, col="red", pch=19, cex=0.7); points(coordsMDS_mod1b, col="blue", pch=19, cex=0.7)
					# points(coordsMDS_mod1c, col="green", pch=19, cex=0.7); points(coordsMDS_mod1d, col="green3", pch=19, cex=0.7)
					# points(coordsMDS_mod2a, col="yellow", pch=19, cex=0.7); points(coordsMDS_mod2b, col="purple", pch=19, cex=0.7)
					# points(coordsMDS_mod2c, col="brown", pch=19, cex=0.7); points(coordsMDS_mod2d, col="pink", pch=19, cex=0.7)
				}
			minNumberOfNApoints = numberOfNApoints[1]; coordsMDS_mod = coordsMDS_modList[[1]]	
			for (j in 2:length(coordsMDS_modList))
				{
					if (numberOfNApoints[j] < minNumberOfNApoints)
						{
							minNumberOfNApoints = numberOfNApoints[j]
							coordsMDS_mod = coordsMDS_modList[[j]]
						}
				}
			dist1 = dist
			dist2 = dist; dist2[,] = 0
			for (j in 2:dim(dist2)[1])
				{
					for (k in 1:j)
						{
							dist2[j,k] = sqrt(((coordsMDS_mod[j,1]-coordsMDS_mod[k,1])^2)+((coordsMDS_mod[j,2]-coordsMDS_mod[k,2])^2))
							dist2[k,j] = dist2[j,k]
						}
				}
			# cor(as.vector(dist1[upper.tri(dist1)]), as.vector(dist2[upper.tri(dist2)]), use="everything", method="pearson")
			# plot(pcoa(dist1)$vectors[,1:2]); dev.new(); plot(pcoa(dist2)$vectors[,1:2], col="red")
			if (showingPlots == TRUE)
				{
					if (i == 1)
						{
							cols = "gray90"
						}	else	{
							cols = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
						}
					plotRaster(envVariables[[i]], cols)
					lines(coordinates[hull,], col="blue", lwd=0.5)
					points(coordinates, cex=0.75)
					points(coordinates, col="black", pch=19, cex=0.7)
					points(coordsMDS_mod, cex=0.8, pch=2)
					points(coordsMDS_mod, col="red", pch=17, cex=0.8)
				}
			tab = cbind(buffer3, coordsMDS_mod)
			if (resistances[[i]] == TRUE) resistance = "R"
			if (resistances[[i]] == FALSE) resistance = "C"
			if (textFile == TRUE)
				{
					sink(file=paste(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),".txt",sep=""))
					for (j in 1:dim(input)[1])
						{
							cat(paste(tab[j,4],tab[j,3],sep="	")); cat("\\n")
						}
					sink(NULL)
				}
			if (fastaAlignment == TRUE)
				{
					fasta = input0
					if (length(row.names(fasta)) > 0)
						{
							fasta = fasta[row.names(fasta)%in%sequences_to_keep,]
							for (j in 1:length(row.names(fasta)))
								{
									row.names(fasta)[j] = paste(unlist(strsplit(tab[j,1],"_"))[1],tab[j,2],tab[j,4],tab[j,3],sep="_")
								}
						}
					if (length(names(fasta)) > 0)
						{
							fasta = fasta[names(fasta)%in%sequences_to_keep,]
							for (j in 1:length(names(fasta)))
								{
									names(fasta)[j] = paste(unlist(strsplit(tab[j,1],"_"))[1],tab[j,2],tab[j,4],tab[j,3],sep="_")
								}
						}	
					# write.table(tab, paste(paste(outputName,"_",rasterNames[i],"_",method,"_",resistance,sep=""),".txt",sep=""), row.names=F, col.names=F, quote=F)
					# write.dna(fasta, paste(paste(outputName,"_",rasterNames[i],"_",method,"_",resistance,sep=""),".fasta",sep=""), format="fasta")
					sink(file=paste(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),".fasta",sep=""))
					for (j in 1:length(row.names(fasta)))
						{
							cat(paste(">",row.names(fasta)[j],sep="")); cat("\\n")
							cat(sequences[j]); cat("\\n")
						}
					sink(NULL)
				}
			if (treeFile == TRUE)
				{
					tips_to_remove = tree$tip.label[!tree$tip.label%in%sequences_to_keep]
					new_tree = drop.tip(tree, tips_to_remove, trim.internal=TRUE, subtree=FALSE, rooted=is.rooted(tree))
					new.labels = c()
					for (j in 1:dim(tab)[1])
						{
							new.labels = c(new.labels, paste(tab[j,1],tab[j,2],tab[j,4],tab[j,3],sep="_"))
						}
					new_tree$tip.label = new.labels
					write.tree(new_tree, paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".tree",sep=""))
					newickTree = scan(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".tree",sep=""), what="", sep="\\n", quiet=TRUE)
					xmlTemplate = scan("Xml_template.xml", what="", sep="\\n", quiet=TRUE)
					xmlBlock = c()
					for (j in 1:dim(tab)[1])
						{
							sequenceName = paste(tab[j,1],tab[j,2],tab[j,4],tab[j,3],sep="_")
							xmlBlock = c(xmlBlock, paste("\	<taxon id=\\"",sequenceName,"\\">",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<date value=\\"",tab[j,2],"\\" direction=\\"forwards\\" units=\\"years\\" precision=\\"0.1\\"/>",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<attr name=\\"lat\\">",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	\	",tab[j,4],sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	</attr>",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<attr name=\\"lon\\">",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	\	",tab[j,3],sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	</attr>",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<!-- START Multivariate diffusion model -->",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<attr name=\\"location\\">",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	\	",tab[j,4]," ",tab[j,3],sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	</attr>",sep=""))
							xmlBlock = c(xmlBlock, paste("\	\	<!-- END Multivariate diffusion model -->",sep=""))
							xmlBlock = c(xmlBlock, paste("\	</taxon>",sep=""))
						}
					sink(file=paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".xml",sep=""))
					cat("<beast>"); cat("\\n")
					cat("\	<taxa id=\\"taxa\\">"); cat("\\n")
					for (j in 1:length(xmlBlock))
						{
							cat(xmlBlock[j]); cat("\\n")
						}
					cat("\	</taxa>"); cat("\\n")
					cat("\\n")
					cat("\	<newick id=\\"startingTree.sim\\">"); cat("\\n")
					cat(newickTree); cat("\\n")
					cat("\	</newick>"); cat("\\n")
					cat("\\n")
					xmlTemplate = gsub("TEMPLATE_NAME",paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),xmlTemplate)
					for (j in 1:length(xmlTemplate))
						{
							cat(xmlTemplate[j]); cat("\\n")
						}
					cat("\\n")
					cat("<beast>"); cat("\\n")
					sink(NULL)
				}
		}
}
