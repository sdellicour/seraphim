treeExtractions <-
function(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores=1) {

\tregisterDoMC(cores=nberOfCores)
\tdir.create(file.path(localTreesDirectory), showWarnings=F)
\tMCCtree = FALSE
\tnewTreeList = list(); count = 1
\tfor (i in 1:length(allTrees)) # to separate trees
\t\t{\t
\t\t\tif (substring(allTrees[[i]],1,4) == "tree")
\t\t\t\t{
\t\t\t\t\tnewTreeList[[count]] = allTrees[[i]]
\t\t\t\t\tcount = count+1
\t\t\t\t}
\t\t}
\tnberOfTrees\t= length(newTreeList)\t
\tnewTreeListMinusBurnIn = newTreeList[(burnIn+1):nberOfTrees]
\tregularSampledTrees = newTreeList[nberOfTreesToSample]
\tif (randomSampling == TRUE)
\t\t{
\t\t\tsampledNewTreeListMinusBurnIn = sample(newTreeListMinusBurnIn,nberOfTreesToSample)
\t\t}\telse\t{
\t\t\tintervale = floor(length(newTreeListMinusBurnIn)/nberOfTreesToSample)
\t\t\tindex = 0
\t\t\tfor (i in 1:nberOfTreesToSample)
\t\t\t\t{
\t\t\t\t\tindex = index + intervale
\t\t\t\t\tregularSampledTrees[i] = newTreeListMinusBurnIn[index]
\t\t\t\t}
\t\t\tsampledNewTreeListMinusBurnIn = regularSampledTrees
\t\t}
\t
\tbuffer = list()
\tbuffer = foreach(i = 1:length(sampledNewTreeListMinusBurnIn)) %dopar% {
\t# for (i in 1:length(sampledNewTreeListMinusBurnIn)) {\t
\t\t\textractionsID = runif(1,0,1)
\t\t\tproperTreeName = paste("ProperTree_", extractionsID, ".tree", sep="")
\t\t\ttree = sampledNewTreeListMinusBurnIn[i]
\t\t\tif (MCCtree == FALSE)
\t\t\t\t{
\t\t\t\t\ttreeID = gsub(pattern="\\\\[.+\\\\];", replacement="", tree)
\t\t\t\t\ttreeID = gsub(pattern=" ", replacement="", treeID)
\t\t\t\t\ttreeID = unlist(strsplit(treeID, "_"))[2]
\t\t\t\t\tcat("Exctracting information from sampled tree n°", treeID, "\\n", sep="");
\t\t\t\t}\telse\t\t{
\t\t\t\t\ttreeID = "MCC_tree"
\t\t\t\t}
\t\t\t# (1) Creation of "tab", a data frame containing coordinates of each node:
\t\t\ttab = gsub(pattern="tree.+\\\\[&R\\\\]", replacement="", tree)\t
\t\t\ttab = gsub(pattern="\\\\}\\\\]:\\\\[&rate=", replacement="\\\\},rate=", tab) 
\t\t\ttab = gsub(pattern="\\\\}\\\\]:\\\\[&+[[:alnum:]]+[[:punct:]]+rate=", replacement="\\\\},rate=", tab)
\t\t\ttab = gsub("\\\\[&CO", "\\\\&&CO", tab)
\t\t\ttab = unlist(strsplit(tab, "\\\\["))[-1]
\t\t\ttab = gsub("&", "", tab)
\t\t\ttab = gsub("\\\\}.+$", "", tab)\t\t
\t\t\ttab = gsub("\\\\{", "", tab)
\t\t\ttab = gsub(paste(".*",coordinateAttributeName,"=",sep=""), "location1=", tab)
\t\t\ttab = gsub(",", ",location2=", tab)
\t    \tfoo = function(x) {x = unlist(strsplit(x, ",")); x}
\t    \ttab = lapply(tab, foo)
\t    \ttab_transit = data.frame()
\t    \tcolNames = unique(gsub("=.+$", "", unlist(tab)))
\t    \tcolNames = gsub(" ", "", colNames)
\t    \tfor (j in seq(tab))
\t    \t\t{
\t    \t\t\tfor (k in seq(colNames))
\t    \t\t\t\t{\t
\t    \t\t\t\t\tind = grep(paste("^", colNames[k], "=", sep=""),tab[[j]])
\t    \t\t\t\t\tif (length(ind) > 0)
\t    \t\t\t\t\t\t{
\t    \t\t\t\t\t\t\tv = as.numeric(gsub(paste(colNames[k], "=", sep=""), "", tab[[j]][ind]))
\t    \t\t\t\t\t\t\ttab_transit[j, k] = v
\t    \t\t\t\t\t\t}
\t    \t\t\t\t}
\t    \t\t}
\t    \ttab = tab_transit
\t    \tcolnames(tab) = colNames   
\t    \tall = which(!is.na(tab[,1]))
\t    \t 
\t\t\t# (2) Removing of all the spatial & mutation rate information:
\t\t\ttree = gsub("\\\\[[^]]*\\\\]", "", tree) # remove all the "[...]"
\t\t\ttree = gsub("tree STATE_.+[[:space:]]+=[[:space:]]+", "", tree)
\t\t\twrite(tree, file=properTreeName)
\t\t\t
\t\t\t# (3) Creation of "mat", a matrix containing coordinates of each node:
\t\t\t\t# tree1 will contain node names (if tip nodes) and branch lengths
\t\t\ttree1 = unlist(strsplit(tree, "\\\\("))
\t   \t\ttree1 = unlist(strsplit(tree1, ","))
\t   \t\ttree1 = unlist(strsplit(tree1, "\\\\)"))
\t\t   \t\t# tree2 will only contain all the branch lengths
\t\t\ttree2 = unlist(strsplit(tree, ":"))[-1]
\t\t\ttree2 = gsub("[( | ) | ;]", "", tree2)\t
\t\t\ttree2 = strsplit(tree2, ",")
\t\t\tfoo = function(x) x = head(x, 1)
\t\t\ttree2 = unlist(lapply(tree2, foo))
\t\t\ttree2 = paste("", tree2, sep=":")
\t\t\ttree2 = c(tree2, ";")
\t\t\tallStats = vector(mode="list", length=dim(tab)[2])
\t\t\tfor (j in seq(allStats))
\t\t\t\t{
\t\t\t\t\ttree3 = tree
\t\t\t\t\tval = tab[,j]
\t\t\t\t\ttree4 = paste(val, tree2, sep="")
\t\t\t\t\t\t# tree4 will only contain one coordinate and the branch length
\t\t\t\t\tfor (k in all)
\t\t\t\t\t\t{
\t\t\t\t\t\t\t# tree3 = gsub(tree1[k], tree4[k], tree3) # PROBLEM
\t\t\t\t\t\t\ttree3 = gsub(paste("\\\\(",tree1[k],sep=""), paste("\\\\(",tree4[k],sep=""), tree3)
\t\t\t\t\t\t\ttree3 = gsub(paste("\\\\)",tree1[k],sep=""), paste("\\\\)",tree4[k],sep=""), tree3)
\t\t\t\t\t\t\ttree3 = gsub(paste(",",tree1[k],sep=""), paste(",",tree4[k],sep=""), tree3)
\t\t\t\t\t\t\t\t# node/tip names are replaced by the coordinate
\t \t\t\t\t\t\ttree5 = read.tree(text = tree3)
\t\t\t\t\t\t\ttree6 = c(tree5$tip.label, tree5$node.label)
\t\t\t\t\t\t\ttree6[tree6 == "NA"] = 9999
\t \t\t\t\t\t\ttree6 = as.numeric(tree6)
\t\t\t\t\t\t\ttree6[tree6 == 9999] = NA
\t\t\t\t\t\t\tallStats[[j]] = tree6
\t\t\t\t\t\t\tnames(allStats)[j] = colnames(tab)[j]\t
\t\t\t\t\t\t}\t# allStats will contain of list of the two coordinates
\t\t\t\t}\t\t\t# but in the correct order (i.e. the "tree$edge order")
\t\t\tmat = matrix(unlist(allStats), ncol=length(allStats), byrow=F)
\t
\t\t\t# (4) Adding the spatial information to the matrix:
\t\t\ttree = read.tree(properTreeName)
\t\t\ttreeEdges = tree$edge
\t\t\ttreeEdgesLengths = tree$edge.length
\t\t\ttreeEdgesLengths = as.matrix(treeEdgesLengths,nrows=(b[1]),ncol=1)
\t\t\ttreeEdgesLengths = cbind(treeEdges,treeEdgesLengths)
\t\t\tmat = cbind(mat,c(1:length(mat[,1])))  
\t\t\tmat = cbind(mat,c(1:length(mat[,1]))) 
\t\t\tcolnames(mat) = c("latitude","longitude","node1","node2")
\t\t\tcolnames(treeEdgesLengths) = c("node1","node2","length")
\t\t\tmat_transit = merge(treeEdgesLengths,mat,"node2")
\t\t\tcolnames(mat_transit) = c("node2","node1","length","endStates","endRate","x")
\t\t\tmat_transit = merge(mat_transit,mat,"node1")
\t\t\tcolnames(mat_transit) = c("node1","node2","length","endLat","endLon","x","startLat","startLon","y")
\t\t\tmat_transit$x = NULL; mat_transit$y = NULL
\t\t\tmat = mat_transit[,c(1,2,3,6,7,4,5)] # changes the order of the columns
\t
\t\t\t# (5) Adding the temporal information to the matrix:
\t\t\tl = length(mat[,1])
\t\t\tif (l > 0)
\t\t\t\t{
\t\t\t\t\tll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
\t\t\t\t\tfor (j in 1:l)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsubMat = mat[j,2]
\t\t\t\t\t\t\tsubMat = subset(mat,mat[,2]==subMat)
\t\t\t\t\t\t\tll[j,1] = subMat[,3]
\t\t\t\t\t\t\tsubMat = subMat[1,1]
\t\t\t\t\t\t\tsubMat1 = subset(mat,mat[,2]==subMat)
\t\t\t\t\t\t\tfor (k in 1:l)
\t\t\t\t\t\t\t\t{
\t\t \t\t\t\t\t\t\tif (nrow(subMat1) > 0)
\t\t \t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tll[j,k+1] = subMat1[,3]
\t \t\t\t\t\t\t\t\t\t\tsubMat2 = subMat1[1,1]
\t \t\t\t\t\t\t\t\t\t\tsubMat1 = subset(mat,mat[,2]==subMat2)
\t \t\t\t\t\t\t\t\t\t}
\t \t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tendNodeL = rowSums(ll)
\t\t \t\t\tmat = cbind(mat, endNodeL) # mat is now the full matrix with the root to node distance for each node
\t\t\t\t\tstartNodeL = matrix(1:l,nrow=l,ncol=1)
\t\t\t\t\tstartNodeL[] = 0
\t\t\t\t\tfor (j in 1:l)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tr = mat[j,1]
\t\t\t\t\t\t\ts = subset(mat,mat[,2]==r)
\t\t\t\t\t\t\tfor (k in 1:l)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (nrow(s) > 0)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tstartNodeL[j,1] = s[,8]
\t \t\t\t\t\t\t\t\t\t}
\t \t\t\t\t\t\t\t}\t
\t\t \t\t\t\t}
\t\t\t\t\tmat = cbind(mat,startNodeL)
\t\t\t\t\tmaxEndLIndice = which.max(mat[,"endNodeL"])
\t \t\t\t\tmaxEndL = mat[maxEndLIndice,"endNodeL"]
\t \t\t\t\tendYear = matrix(mat[,"endNodeL"]-maxEndL)
\t \t\t\t\tendYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
\t\t\t\t\tstartYear = matrix(mat[,"startNodeL"]-maxEndL)
\t \t\t\t\tstartYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
\t \t\t\t\tmat = cbind(mat,startYear,endYear)
\t \t\t\t\t
\t \t\t\t\t# (6) Adding surface distance and tree ID:
\t \t\t\t\t# distances = matrix(0, nrow=dim(mat)[1], ncol=1)
\t \t\t\t\t# colnames(distances) = c("greatCircleDist_km")
\t \t\t\t\ttreeIDs = matrix(0, nrow=dim(mat)[1], ncol=1)
\t \t\t\t\tcolnames(treeIDs) = c("treeID")
\t \t\t\t\tfor (j in 1:dim(mat)[1])
\t \t\t\t\t\t{
\t \t\t\t\t\t\t# x1 = cbind(mat[j,"startLon"],mat[j,"startLat"])
\t \t\t\t\t\t\t# x2 = cbind(mat[j,"endLon"],mat[j,"endLat"])
\t \t\t\t\t\t\t# distances[j] = rdist.earth(x1, x2, miles=F, R=NULL)
\t \t\t\t\t\t\ttreeIDs[j] = treeID
\t \t\t\t\t\t}
\t \t\t\t\t# mat = cbind(mat,distances,treeIDs)
\t \t\t\t\tmat = cbind(mat,treeIDs)
\t \t\t
\t \t\t\t\tfile.remove(properTreeName)
\t\t\t\t\tfile = as.character(paste(localTreesDirectory, "/TreeExtractions_", i, ".csv", sep=""))
\t\t\t\t\twrite.csv(mat, file, row.names=F, quote=F)
\t\t\t\t}\telse\t{
\t\t\t\t\tfile.remove(properTreeName)
\t\t\t\t\tprint(paste("A problem occured with tree extraction ",i,sep=""))
\t\t\t\t}
\t\t}
}
