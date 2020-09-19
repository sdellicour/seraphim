simulatorRRW1 <-
function(tree, rates, sigmas=c(0.1,0.1), cor=0, envVariables=list(), mostRecentSamplingDatum,
\t\t\t\t\t\t ancestPosition, reciprocalRates=TRUE, n1=100, n2=0, showingPlots=FALSE, newPlot=TRUE, fixedNodes=c()) {
\t
\trotation = function(pt1, pt2, angle)
\t\t{
\t\t\ts = sin(angle); c = cos(angle)
\t\t\tx = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
\t\t\tx_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
\t\t\tx_new = x_new+pt1[1]; y_new = y_new+pt1[2]
\t\t\treturn(c(x_new,y_new))
\t\t}
\tnodesOnly = FALSE; pointCol = "red"
\tcolNames= c("node1","node2","length","startLat","startLon","endLat","endLon",
\t\t\t\t"endNodeL","startNodeL","startYear","endYear","greatCircleDist_km")
\tsimulation = matrix(nrow=length(tree$edge.length), ncol=length(colNames))
\tcolnames(simulation) = colNames
\t# if (model == "fixed") phi_b = rep(1, length(tree$edge.length))
\t# if (model == "cauchy") phi_b = rgamma(length(tree$edge.length), shape=0.5, scale=0.5)
\t# if (model == "gamma") phi_b = rgamma(length(tree$edge.length), shape=halfDF, scale=halfDF)
\t# if (model == "logN") phi_b = rlnorm(length(tree$edge.length), meanlog=1, sdlog=sdLogN)
\t# if (model == "cauchy") sd_BM = sqrt(tree$edge.length/phi_b) # corresponds to BEAST reciprocalRates="true"
\t# if (model == "gamma") sd_BM = sqrt(tree$edge.length/phi_b) # corresponds to BEAST reciprocalRates="true"
\t# if (model == "logN") sd_BM = sqrt(tree$edge.length*phi_b) # corresponds to BEAST reciprocalRates="false"
\tphi_b = rates
\tif (reciprocalRates == FALSE) sd_BM = sqrt(tree$edge.length*phi_b)
\tif (reciprocalRates == TRUE) sd_BM = sqrt(tree$edge.length/phi_b)
\tnd = node.depth(tree); nd_max = max(nd)
\tt1 = rep(ancestPosition[2], length(tree$tip.label)+tree$Nnode)
\tt2 = rep(ancestPosition[1], length(tree$tip.label)+tree$Nnode)
\tif (showingPlots == TRUE)
\t\t{
\t\t\tif ((newPlot == TRUE)&(length(envVariables) > 0))
\t\t\t\t{
\t\t\t\t\tplotRaster(rast=envVariables[[1]], cols="gray90", colNA="white", addBox=F)
\t\t\t\t}
\t\t\tpoints(cbind(ancestPosition[1],ancestPosition[2]), pch=16, col=pointCol, cex=0.5)
\t\t}
\ti = 0
\twhile (i != (nd_max-1))
\t\t{
\t\t\ti = i+1
\t\t\tmy_nodes = which(nd==nd_max-i)
\t\t\tif (length(my_nodes) > 0)
\t\t\t\t{
\t\t\t\t\tfor (j in 1:length(my_nodes))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tmy_node = my_nodes[j]; # print(c(i,my_nodes[j]))
\t\t\t\t\t\t\tparent_branch = match(my_node, tree$edge[,2])
\t\t\t\t\t\t\tparent_node = tree$edge[parent_branch,1] 
\t\t\t\t\t\t\tsimulatingNode = TRUE
\t\t\t\t\t\t\tif (my_node%in%fixedNodes)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tsimulatingNode = FALSE
\t\t\t\t\t\t\t\t\tindex = which(tree$edge[,2] == my_node)
\t\t\t\t\t\t\t\t\tnew_t1 = tree$annotations[[index]]$location[[1]]
\t\t\t\t\t\t\t\t\tnew_t2 = tree$annotations[[index]]$location[[2]]
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tif (simulatingNode == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tonTheArea = FALSE
\t\t\t\t\t\t\t\t\tsd_bm = sd_BM[parent_branch]
\t\t\t\t\t\t\t\t\tincrement1 = rnorm(1)*sd_bm
\t\t\t\t\t\t\t\t\tincrement2 = rnorm(1)*sd_bm
\t\t\t\t\t\t\t\t\tincrement2 = sigmas[2]*((cor*increment1)+(sqrt(1-cor^2)*increment2))
\t\t\t\t\t\t\t\t\tincrement1 = sigmas[1]*increment1
\t\t\t\t\t\t\t\t\tnew_t1 = t1[parent_node] + increment1 
\t\t\t\t\t\t\t\t\tnew_t2 = t2[parent_node] + increment2
\t\t\t\t\t\t\t\t\tif (length(envVariables) > 0)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tonTheArea = TRUE
\t\t\t\t\t\t\t\t\t\t\tfor (k in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(envVariables[[k]],cbind(new_t2,new_t1))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tonTheArea = TRUE
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tif (onTheArea == FALSE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tc2 = 0; c1 = 0
\t\t\t\t\t\t\t\t\t\t\tpt1 = cbind(t2[parent_node],t1[parent_node])
\t\t\t\t\t\t\t\t\t\t\tpt2 = cbind(new_t2,new_t1)
\t\t\t\t\t\t\t\t\t\t\twhile (onTheArea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tc2 = c2+1; # print(c(c2,c1))
\t\t\t\t\t\t\t\t\t\t\t\t\tif (n1 > 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (k in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(envVariables[[k]],cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnew_t1 = pt2_rotated[2]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnew_t2 = pt2_rotated[1]
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tif (c2 > n1)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tc2 = 0; c1 = c1+1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (c1 > n2)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = TRUE; i = 0
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t# print(paste0("...re-simulating a branch - node:", my_node))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tincrement1 = rnorm(1)*sd_bm
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tincrement2 = rnorm(1)*sd_bm
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tincrement2 = sigmas[2]*((cor*increment1)+(sqrt(1-cor^2)*increment2))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tincrement1 = sigmas[1]*increment1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnew_t1 = t1[parent_node] + increment1 
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnew_t2 = t2[parent_node] + increment2
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = TRUE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfor (k in 1:length(envVariables))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(envVariables[[k]],cbind(new_t2,new_t1))))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tonTheArea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tt1[my_node] = new_t1\t
\t\t\t\t\t\t\tt2[my_node] = new_t2
\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (nodesOnly == FALSE)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tsegments(t2[parent_node], t1[parent_node], new_t2, new_t1, col=pointCol, lwd=0.2)
\t\t\t\t\t\t\t\t\t\t\tpoints(cbind(new_t2,new_t1), pch=16, col=pointCol, cex=0.25)
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tpoints(cbind(new_t2,new_t1), pch=16, col=pointCol, cex=0.25)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (i == 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste0("...re-starting the simulation\\n"))
\t\t\t\t\t\t\tt1 = rep(ancestPosition[2], length(tree$tip.label)+tree$Nnode)
\t\t\t\t\t\t\tt2 = rep(ancestPosition[1], length(tree$tip.label)+tree$Nnode)
\t\t\t\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tplotRaster(rast=envVariables[[1]], cols="gray90", colNA="white", addBox=F, new=F)
\t\t\t\t\t\t\t\t\tpoints(cbind(ancestPosition[1],ancestPosition[2]), pch=16, col=pointCol, cex=0.5)
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t}
\tx = t2; y = t1
\tfor (i in 1:dim(tree$edge)[1])
\t\t{
\t\t\tnode_i = tree$edge[i,1]
\t\t\tnode_f = tree$edge[i,2]
\t\t\tsimulation[i,"node1"] = node_i
\t\t\tsimulation[i,"node2"] = node_f
\t\t\tsimulation[i,"length"] = tree$edge.length[i]
\t\t\tsimulation[i,"startLat"] = y[node_i]
\t\t\tsimulation[i,"startLon"] = x[node_i]
\t\t\tsimulation[i,"endLat"] = y[node_f]
\t\t\tsimulation[i,"endLon"] = x[node_f]
\t\t\tx1 = cbind(x[node_i],y[node_i]); x2 = cbind(x[node_f],y[node_f])
\t\t\tsimulation[i,"greatCircleDist_km"] = rdist.earth(x1, x2, miles=FALSE, R=NULL)
\t\t}\t
\tl = length(simulation[,1])
\tll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
\tfor (j in 1:l)
\t\t{
\t\t\tsubMat = simulation[j,2]
\t\t\tsubMat = subset(simulation,simulation[,2]==subMat)
\t\t\tll[j,1] = subMat[,3]
\t\t\tsubMat = subMat[1,1]
\t\t\tsubMat1 = subset(simulation,simulation[,2]==subMat)
\t\t\tfor (k in 1:l)
\t\t\t\t{
\t \t\t\t\tif (nrow(subMat1) > 0)
\t \t\t\t\t\t{
\t\t\t\t\t\t\tll[j,k+1] = subMat1[,3]
 \t\t\t\t\t\t\tsubMat2 = subMat1[1,1]
 \t\t\t\t\t\t\tsubMat1 = subset(simulation,simulation[,2]==subMat2)
 \t\t\t\t\t\t}
 \t\t\t\t}
\t\t}
\tendNodeL = rowSums(ll) # root to node distance for each node
\tsimulation[,"endNodeL"] = endNodeL
\tstartNodeL = matrix(1:l,nrow=l,ncol=1)
\tstartNodeL[] = 0
\tfor (j in 1:l)
\t\t{
\t\t\tr = simulation[j,1]
\t\t\ts = subset(simulation,simulation[,2]==r)
\t\t\tfor (k in 1:l)
\t\t\t\t{
\t\t\t\t\tif (nrow(s) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tstartNodeL[j,1] = s[,"endNodeL"]
 \t\t\t\t\t\t}
 \t\t\t\t}\t
\t \t}
\tsimulation[,"startNodeL"] = startNodeL
\tmaxEndLIndice = which.max(simulation[,"endNodeL"])
 \tmaxEndL = simulation[maxEndLIndice,"endNodeL"]
 \tendYear = matrix(simulation[,"endNodeL"]-maxEndL)
 \tendYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
\tstartYear = matrix(simulation[,"startNodeL"]-maxEndL)
 \tstartYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
 \tsimulation[,c("startYear","endYear")] = cbind(startYear,endYear)
 \tif (showingPlots == TRUE) dev.off()
\treturn(simulation)
}
