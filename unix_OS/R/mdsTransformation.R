mdsTransformation <-
function(input, envVariables=list(), pathModel=2, resistances=c(T), avgResistances=c(T), fourCells=F, outputName="", OS="Unix")\t{

\tinput0 = input; showingPlots = FALSE
\tnberOfCores = 1; registerDoMC(cores=nberOfCores)
\ttextFile = FALSE; fastaAlignment = FALSE; treeFile = FALSE
\tif (attr(input,"class") == "data.frame") textFile = TRUE
\tif (attr(input,"class") == "DNAbin") fastaAlignment = TRUE
\tif (attr(input,"class") == "phylo") treeFile = TRUE
\tif (textFile == TRUE)
\t\t{
\t\t\tsequencesID = seq(1,dim(input)[1],1); years = rep(0,dim(input)[1]); locations1 = input[,2]; locations2 = input[,1]
\t\t\tcoordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
\t\t}
\tif (fastaAlignment == TRUE)
\t\t{
\t\t\tfastaChar = as.character(input)
\t\t\tlabels = rownames(fastaChar)
\t\t\tnames = matrix(nrow=length(labels), ncol=4)
\t\t\tsequences = matrix(nrow=dim(fastaChar)[1], ncol=1)
\t\t\tfor (k in 1:dim(fastaChar)[1])
\t\t\t\t{
\t\t\t\t\tnames[k,1:4] = unlist(strsplit(labels[k], "_"))
\t\t\t\t\tsequence = ""
\t\t\t\t\tfor (l in 1:dim(fastaChar)[2]) sequence = paste(sequence,fastaChar[k,l],sep="")
\t\t\t\t\tsequence = gsub("a","A",sequence); sequence = gsub("c","C",sequence)
\t\t\t\t\tsequence = gsub("g","G",sequence); sequence = gsub("t","T",sequence)
\t\t\t\t\tsequences[k,1] = sequence
\t\t\t\t}
\t\t\tsequencesID = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
\t\t\tcoordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
\t\t}
\tif (treeFile == TRUE)
\t\t{
\t\t\ttree = input
\t\t\tlabels = tree$tip.label
\t\t\tlabels = gsub("'","",labels)
\t\t\tnames = matrix(nrow=length(labels), ncol=4)
\t\t\tfor (k in 1:length(labels))
\t\t\t\t{
\t\t\t\t\tnames[k,1:4] = unlist(strsplit(labels[k], "_"))
\t\t\t\t}
\t\t\tsequencesID = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
\t\t\tcoordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
\t\t}
\tnullRaster = envVariables[[1]]; nullRaster[!is.na(nullRaster[])] = 1
\tnames(nullRaster) = "nullRaster"
\tbuffer1 = list(nullRaster)
\tbuffer2 = list(T); buffer3 = list(T)
\tfor (i in 1:length(envVariables))
\t\t{
\t\t\tbuffer1[[i+1]] = envVariables[[i]]
\t\t\tbuffer2[[i+1]] = resistances[[i]]
\t\t\tbuffer3[[i+1]] = avgResistances[[i]]
\t\t}
\tenvVariables = list(); envVariables = buffer1
\tresistances = list(); resistances = buffer2
\tavgResistances = list(); avgResistances = buffer3
\trasterNames = c()
\tfor (i in 1:length(envVariables)) rasterNames = c(rasterNames, names(envVariables[[i]]))
\tif (pathModel == 3)
\t\t{
\t\t\tif("CS_envVariables"%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), "CS_envVariables"))
\t\t\tfor (i in 1:length(envVariables))
\t\t\t\t{
\t\t\t\t\t# name = paste("CS_envVariables/", gsub("_0.","_0,",names(envVariables[[i]])), "_cs.asc", sep="")
\t\t\t\t\tname = paste("CS_envVariables/", names(envVariables[[i]]), "_cs.asc", sep="")
\t\t\t\t\twriteRaster(envVariables[[i]], name, overwrite=T)
\t\t\t\t}
\t\t}
\tdists = list(); euclDis = matrix(0, nrow=dim(coordinates)[1], ncol=dim(coordinates)[1])
\tbuffer1 = coordinates[,1:2]; coordinates0 = coordinates
\tcoordinates = cbind(as.numeric(coordinates0[,3]),as.numeric(coordinates0[,4]))
\tminEuclDis = sqrt(((coordinates[2,1]-coordinates[1,1])^2)+((coordinates[2,2]-coordinates[1,2])^2))
\tmaxEuclDis = sqrt(((coordinates[2,1]-coordinates[1,1])^2)+((coordinates[2,2]-coordinates[1,2])^2))
\tfor (j in 2:dim(coordinates)[1])
\t\t{
\t\t\tfor (k in 1:j)
\t\t\t\t{
\t\t\t\t\td = sqrt(((coordinates[j,1]-coordinates[k,1])^2)+((coordinates[j,2]-coordinates[k,2])^2))
\t\t\t\t\teuclDis[k,j] = d; euclDis[j,k] = d
\t\t\t\t\tif (minEuclDis > d) minEuclDis = d
\t\t\t\t\tif (maxEuclDis < d) maxEuclDis = d
\t\t\t\t}
\t\t}
\thullRaster = nullRaster
\thull = chull(coordinates)
\thull = c(hull,hull[1])
\tp = Polygon(coordinates[hull,])
\tps = Polygons(list(p),1)
\tsps = SpatialPolygons(list(ps))
\tpointsRaster = rasterize(coordinates, crop(hullRaster, sps, snap="out"))
\tpointsRaster[!is.na(pointsRaster[])] = 0
\thullRaster = crop(hullRaster, sps, snap="out")
\tbufferRaster = hullRaster
\thullRaster = mask(hullRaster, sps, snap="out")
\thullRaster[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
\tbuffer = list()
\tbuffer = foreach(i = 1:length(envVariables)) %dopar% {
\t# for (i in 1:length(envVariables)) {
\t\t\tsequences_to_keep = labels
\t\t\tdist = matrix(0, nrow=length(coordinates[,1]), ncol=length(coordinates[,1]))
\t\t\tif (pathModel == 1)
\t\t\t\t{
\t\t\t\t\tmethod = "SL"
\t\t\t\t\tlinesList = list()
\t\t\t\t\tc = 0
\t\t\t\t\tfor (j in 2:length(coordinates[,1]))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfor (k in 1:(j-1))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\tpoints = rbind(coordinates[k,], coordinates[j,])
\t\t\t\t\t\t\t\t\tlinesList[[c]] = Lines(list(Line(points)),c)
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tlines = SpatialLines(linesList)
\t\t\t\t\textractions = extract(envVariables[[i]], lines)
\t\t\t\t\tc = 0
\t\t\t\t\tfor (j in 2:length(coordinates[,1]))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfor (k in 1:(j-1))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tc = c+1
\t\t\t\t\t\t\t\t\tdist[k,j] = sum(extractions[[c]], na.rm=T)
\t\t\t\t\t\t\t\t\tdist[j,k] = dist[k,j]
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tif (pathModel == 2)
\t\t\t\t{
\t\t\t\t\tmethod = "LC"
\t\t\t\t\tif (resistances[[i]] == T)
\t\t\t\t\t\t{
\t\t\t\t\t\t\ttrEnvVariable = transition(envVariables[[i]], function(x) 1/mean(x), directions=4)
\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\ttrEnvVariable = transition(envVariables[[i]], mean, directions=4)
\t\t\t\t\t\t}
\t\t\t\t\ttrEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T)
\t\t\t\t\tdist = costDistance(trEnvVariableCorr, as.matrix(coordinates), as.matrix(coordinates))
\t\t\t\t}
\t\t\tif (pathModel == 3)
\t\t\t\t{
\t\t\t\t\tmethod = "RW"\t\t\t\t
\t\t\t\t\tdist = circuitScape1(envVariables[[i]], paste("CS_envVariables/",names(envVariables[[i]]),"_cs",sep=""), 
\t\t\t\t\t\t\t\t\t\tresistances[[i]], avgResistances[[i]], fourCells, coordinates, coordinates, OS, outputName, i)
\t\t\t\t\tdist = as.matrix(dist)
\t\t\t\t}
\t\t\tdist[!is.finite(dist[])] = NA
\t\t\t# To remove lines/columns of NA:
\t\t\tlines_to_remove = c()
\t\t\tfor (j in 1:dim(dist)[1])
\t\t\t\t{
\t\t\t\t\tallNA = TRUE
\t\t\t\t\tnas = 0
\t\t\t\t\tfor (k in 1:dim(dist)[2])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (!is.na(dist[j,k])) 
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (dist[j,k]!=0) allNA = FALSE
\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\tnas = nas+1
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (nas >= 3) allNA = TRUE
\t\t\t\t\tif (allNA == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tlines_to_remove = c(lines_to_remove, j)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tif (length(lines_to_remove) > 0)
\t\t\t\t{
\t\t\t\t\tcat(paste("Removed sequences:","\\n",sep=""))
\t\t\t\t\tfor (j in 1:length(lines_to_remove))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste("   ", buffer1[j,1],"\\n",sep=""))\t
\t\t\t\t\t\t}
\t\t\t\t\tdist = dist[-lines_to_remove,]; dist = dist[,-lines_to_remove]
\t\t\t\t\tbuffer2 = buffer1[-lines_to_remove,]
\t\t\t\t\tsequences_to_keep = sequences_to_keep[-lines_to_remove]
\t\t\t\t}\telse\t{
\t\t\t\t\tbuffer2 = buffer1
\t\t\t\t}
\t\t\t# To remove lines/columns with a few or a single NA value(s):
\t\t\tlines_to_remove = c(); nberOfNA = c()
\t\t\tfor (j in 1:dim(dist)[1])
\t\t\t\t{
\t\t\t\t\tatLeastOneNA = FALSE; n = 0
\t\t\t\t\tcolIDs = c()
\t\t\t\t\tfor (k in 1:dim(dist)[2])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif (is.na(dist[j,k]))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tatLeastOneNA = TRUE
\t\t\t\t\t\t\t\t\tn = n+1
\t\t\t\t\t\t\t\t\tcolIDs = c(colIDs, k)
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\t# if (length(colIDs) > 0) print(colIDs)
\t\t\t\t\tif (atLeastOneNA == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tlines_to_remove = c(lines_to_remove, j)
\t\t\t\t\t\t\tnberOfNA = c(nberOfNA, n)
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tif (length(lines_to_remove) > 0)
\t\t\t\t{
\t\t\t\t\tlines_to_remove = lines_to_remove[order(nberOfNA,decreasing=T)]
\t\t\t\t\tstillNA = TRUE; J = 0
\t\t\t\t\twhile (stillNA == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tJ = J+1
\t\t\t\t\t\t\ttest = dist[-lines_to_remove[1:J],-lines_to_remove[1:J]]
\t\t\t\t\t\t\tatLeastOneNA = FALSE; n = 0
\t\t\t\t\t\t\tfor (j in 1:dim(test)[1])
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tfor (k in 1:dim(test)[2])
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (is.na(test[j,k]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tatLeastOneNA = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\tn = n+1
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tprint(n)\t
\t\t\t\t\t\t\tif (atLeastOneNA == FALSE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tdist = dist[-lines_to_remove[1:J],]; dist = dist[,-lines_to_remove[1:J]]
\t\t\t\t\t\t\t\t\tbuffer3 = buffer2[-lines_to_remove[1:J],]
\t\t\t\t\t\t\t\t\tsequences_to_keep = sequences_to_keep[-lines_to_remove]
\t\t\t\t\t\t\t\t\tcat(paste("Removed sequences:","\\n",sep=""))
\t\t\t\t\t\t\t\t\tfor (j in 1:length(lines_to_remove[1:J]))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tcat(paste("   ", buffer1[lines_to_remove[j],1],"\\n",sep=""))\t
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tstillNA = FALSE
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}\telse\t{
\t\t\t\t\tbuffer3 = buffer2
\t\t\t\t}
\t\t\tminDist = min(dist, na.rm=T); maxDist = max(dist, na.rm=T)
\t\t\tdist_mod = (((dist-minDist)/(maxDist-minDist))*(maxEuclDis-minEuclDis))+minEuclDis
\t\t\tcoordsMDS = cmdscale(dist_mod, eig=T, k=2); coordsMDS_mod = coordsMDS$points
\t\t\tcoordsMDS_mod[,1] = coordsMDS_mod[,1]-median(coordsMDS_mod[,1])
\t\t\tcoordsMDS_mod[,2] = coordsMDS_mod[,2]-median(coordsMDS_mod[,2])
\t\t\tcoordsMDS_mod1 = coordsMDS_mod; coordsMDS_mod2 = coordsMDS_mod
\t\t\tcoordsMDS_mod1[,1] = coordsMDS_mod1[,1]+median(coordinates[,1])
\t\t\tcoordsMDS_mod1[,2] = coordsMDS_mod1[,2]+median(coordinates[,2])
\t\t\tcoordsMDS_mod2 = coordsMDS_mod2[,2:1]
\t\t\tcoordsMDS_mod2[,1] = coordsMDS_mod2[,1]+median(coordinates[,1])
\t\t\tcoordsMDS_mod2[,2] = coordsMDS_mod2[,2]+median(coordinates[,2])
\t\t\tcoordsMDS_mod1a = coordsMDS_mod1; coordsMDS_mod1b = coordsMDS_mod1; coordsMDS_mod1c = coordsMDS_mod1; coordsMDS_mod1d = coordsMDS_mod1
\t\t\tcoordsMDS_mod2a = coordsMDS_mod2; coordsMDS_mod2b = coordsMDS_mod2; coordsMDS_mod2c = coordsMDS_mod2; coordsMDS_mod2d = coordsMDS_mod2
\t\t\tcoordsMDS_mod1b[,1] = coordsMDS_mod1b[,1]+(2*(mean(coordsMDS_mod1b[,1])-coordsMDS_mod1b[,1]))
\t\t\tcoordsMDS_mod1c[,2] = coordsMDS_mod1c[,2]+(2*(mean(coordsMDS_mod1c[,2])-coordsMDS_mod1c[,2]))
\t\t\tcoordsMDS_mod1d[,1] = coordsMDS_mod1d[,1]+(2*(mean(coordsMDS_mod1d[,1])-coordsMDS_mod1d[,1]))
\t\t\tcoordsMDS_mod1d[,2] = coordsMDS_mod1d[,2]+(2*(mean(coordsMDS_mod1d[,2])-coordsMDS_mod1d[,2]))
\t\t\tcoordsMDS_mod2b[,1] = coordsMDS_mod2b[,1]+(2*(mean(coordsMDS_mod2b[,1])-coordsMDS_mod2b[,1]))
\t\t\tcoordsMDS_mod2c[,2] = coordsMDS_mod2c[,2]+(2*(mean(coordsMDS_mod2c[,2])-coordsMDS_mod2c[,2]))
\t\t\tcoordsMDS_mod2d[,1] = coordsMDS_mod2d[,1]+(2*(mean(coordsMDS_mod2d[,1])-coordsMDS_mod2d[,1]))
\t\t\tcoordsMDS_mod2d[,2] = coordsMDS_mod2d[,2]+(2*(mean(coordsMDS_mod2d[,2])-coordsMDS_mod2d[,2]))
\t\t\tcoordsMDS_modList = list(coordsMDS_mod1a,coordsMDS_mod1b,coordsMDS_mod1c,coordsMDS_mod1d,coordsMDS_mod2a,coordsMDS_mod2b,coordsMDS_mod2c,coordsMDS_mod2d)
\t\t\tnumberOfNApoints = c()
\t\t\tfor (j in 1:length(coordsMDS_modList))
\t\t\t\t{
\t\t\t\t\tn = 0
\t\t\t\t\tfor (k in 1:dim(coordsMDS_modList[[j]])[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tif(is.na(raster::extract(hullRaster,cbind(coordsMDS_modList[[j]][k,1],coordsMDS_modList[[j]][k,2])))) n = n+1
\t\t\t\t\t\t}
\t\t\t\t\tnumberOfNApoints = c(numberOfNApoints, n)
\t\t\t\t}
\t\t\tif (showingPlots == TRUE)
\t\t\t\t{\t
\t\t\t\t\t# plotRaser(nullRaster, col="grey90")
\t\t\t\t\t# lines(coordinates[hull,], lwd=0.5); points(coordinates, pch=19, cex=0.7)
\t\t\t\t\t# points(coordsMDS_mod1a, col="red", pch=19, cex=0.7); points(coordsMDS_mod1b, col="blue", pch=19, cex=0.7)
\t\t\t\t\t# points(coordsMDS_mod1c, col="green", pch=19, cex=0.7); points(coordsMDS_mod1d, col="green3", pch=19, cex=0.7)
\t\t\t\t\t# points(coordsMDS_mod2a, col="yellow", pch=19, cex=0.7); points(coordsMDS_mod2b, col="purple", pch=19, cex=0.7)
\t\t\t\t\t# points(coordsMDS_mod2c, col="brown", pch=19, cex=0.7); points(coordsMDS_mod2d, col="pink", pch=19, cex=0.7)
\t\t\t\t}
\t\t\tminNumberOfNApoints = numberOfNApoints[1]; coordsMDS_mod = coordsMDS_modList[[1]]\t
\t\t\tfor (j in 2:length(coordsMDS_modList))
\t\t\t\t{
\t\t\t\t\tif (numberOfNApoints[j] < minNumberOfNApoints)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tminNumberOfNApoints = numberOfNApoints[j]
\t\t\t\t\t\t\tcoordsMDS_mod = coordsMDS_modList[[j]]
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tdist1 = dist
\t\t\tdist2 = dist; dist2[,] = 0
\t\t\tfor (j in 2:dim(dist2)[1])
\t\t\t\t{
\t\t\t\t\tfor (k in 1:j)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tdist2[j,k] = sqrt(((coordsMDS_mod[j,1]-coordsMDS_mod[k,1])^2)+((coordsMDS_mod[j,2]-coordsMDS_mod[k,2])^2))
\t\t\t\t\t\t\tdist2[k,j] = dist2[j,k]
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\t# cor(as.vector(dist1[upper.tri(dist1)]), as.vector(dist2[upper.tri(dist2)]), use="everything", method="pearson")
\t\t\t# plot(pcoa(dist1)$vectors[,1:2]); dev.new(); plot(pcoa(dist2)$vectors[,1:2], col="red")
\t\t\tif (showingPlots == TRUE)
\t\t\t\t{
\t\t\t\t\tif (i == 1)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcols = "gray90"
\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
\t\t\t\t\t\t}
\t\t\t\t\tplotRaster(envVariables[[i]], cols)
\t\t\t\t\tlines(coordinates[hull,], col="blue", lwd=0.5)
\t\t\t\t\tpoints(coordinates, cex=0.75)
\t\t\t\t\tpoints(coordinates, col="black", pch=19, cex=0.7)
\t\t\t\t\tpoints(coordsMDS_mod, cex=0.8, pch=2)
\t\t\t\t\tpoints(coordsMDS_mod, col="red", pch=17, cex=0.8)
\t\t\t\t}
\t\t\ttab = cbind(buffer3, coordsMDS_mod)
\t\t\tif (resistances[[i]] == TRUE) resistance = "R"
\t\t\tif (resistances[[i]] == FALSE) resistance = "C"
\t\t\tif (textFile == TRUE)
\t\t\t\t{
\t\t\t\t\tsink(file=paste(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),".txt",sep=""))
\t\t\t\t\tfor (j in 1:dim(input)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste(tab[j,4],tab[j,3],sep="\t")); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t\tif (fastaAlignment == TRUE)
\t\t\t\t{
\t\t\t\t\tfasta = input0
\t\t\t\t\tif (length(row.names(fasta)) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfasta = fasta[row.names(fasta)%in%sequences_to_keep,]
\t\t\t\t\t\t\tfor (j in 1:length(row.names(fasta)))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\trow.names(fasta)[j] = paste(unlist(strsplit(tab[j,1],"_"))[1],tab[j,2],tab[j,4],tab[j,3],sep="_")
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (length(names(fasta)) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfasta = fasta[names(fasta)%in%sequences_to_keep,]
\t\t\t\t\t\t\tfor (j in 1:length(names(fasta)))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tnames(fasta)[j] = paste(unlist(strsplit(tab[j,1],"_"))[1],tab[j,2],tab[j,4],tab[j,3],sep="_")
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}\t
\t\t\t\t\t# write.table(tab, paste(paste(outputName,"_",rasterNames[i],"_",method,"_",resistance,sep=""),".txt",sep=""), row.names=F, col.names=F, quote=F)
\t\t\t\t\t# write.dna(fasta, paste(paste(outputName,"_",rasterNames[i],"_",method,"_",resistance,sep=""),".fasta",sep=""), format="fasta")
\t\t\t\t\tsink(file=paste(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),".fasta",sep=""))
\t\t\t\t\tfor (j in 1:length(row.names(fasta)))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste(">",row.names(fasta)[j],sep="")); cat("\\n")
\t\t\t\t\t\t\tcat(sequences[j]); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t\tif (treeFile == TRUE)
\t\t\t\t{
\t\t\t\t\ttips_to_remove = tree$tip.label[!tree$tip.label%in%sequences_to_keep]
\t\t\t\t\tnew_tree = drop.tip(tree, tips_to_remove, trim.internal=TRUE, subtree=FALSE, rooted=is.rooted(tree))
\t\t\t\t\tnew.labels = c()
\t\t\t\t\tfor (j in 1:dim(tab)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tnew.labels = c(new.labels, paste(tab[j,1],tab[j,2],tab[j,4],tab[j,3],sep="_"))
\t\t\t\t\t\t}
\t\t\t\t\tnew_tree$tip.label = new.labels
\t\t\t\t\twrite.tree(new_tree, paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".tree",sep=""))
\t\t\t\t\tnewickTree = scan(paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".tree",sep=""), what="", sep="\\n", quiet=TRUE)
\t\t\t\t\txmlTemplate = scan("Xml_template.xml", what="", sep="\\n", quiet=TRUE)
\t\t\t\t\txmlBlock = c()
\t\t\t\t\tfor (j in 1:dim(tab)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsequenceName = paste(tab[j,1],tab[j,2],tab[j,4],tab[j,3],sep="_")
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t<taxon id=\\"",sequenceName,"\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<date value=\\"",tab[j,2],"\\" direction=\\"forwards\\" units=\\"years\\" precision=\\"0.1\\"/>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lat\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",tab[j,4],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lon\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",tab[j,3],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- START Multivariate diffusion model -->",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"location\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",tab[j,4]," ",tab[j,3],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- END Multivariate diffusion model -->",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t</taxon>",sep=""))
\t\t\t\t\t\t}
\t\t\t\t\tsink(file=paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,".xml",sep=""))
\t\t\t\t\tcat("<beast>"); cat("\\n")
\t\t\t\t\tcat("\\t<taxa id=\\"taxa\\">"); cat("\\n")
\t\t\t\t\tfor (j in 1:length(xmlBlock))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(xmlBlock[j]); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tcat("\\t</taxa>"); cat("\\n")
\t\t\t\t\tcat("\\n")
\t\t\t\t\tcat("\\t<newick id=\\"startingTree.sim\\">"); cat("\\n")
\t\t\t\t\tcat(newickTree); cat("\\n")
\t\t\t\t\tcat("\\t</newick>"); cat("\\n")
\t\t\t\t\tcat("\\n")
\t\t\t\t\txmlTemplate = gsub("TEMPLATE_NAME",paste(outputName,"_MDS_",rasterNames[i],"_",method,"_",resistance,sep=""),xmlTemplate)
\t\t\t\t\tfor (j in 1:length(xmlTemplate))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(xmlTemplate[j]); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tcat("\\n")
\t\t\t\t\tcat("<beast>"); cat("\\n")
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t}
}
