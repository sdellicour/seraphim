cartogramTransformation <-
function(input, envVariables=list(), resistances=c(T), outputName="")\t{
\t
\tinput0 = input; showingPlots = FALSE
\ttextFile = FALSE; fastaAlignment = FALSE; treeFile = FALSE
\tfor (i in 1:length(envVariables))
\t\t{
\t\t\tif (resistances[[i]] == FALSE)
\t\t\t\t{
\t\t\t\t\tenvVariables[[i]][!is.na(envVariables[[i]][])] = 1/envVariables[[i]][!is.na(envVariables[[i]][])]
\t\t\t\t}
\t\t}
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
\t\t\tsequencesID = rownames(fastaChar)
\t\t\tnames = matrix(nrow=dim(fastaChar)[1], ncol=4)
\t\t\tsequences = matrix(nrow=dim(fastaChar)[1], ncol=1)
\t\t\tfor (k in 1:dim(fastaChar)[1])
\t\t\t\t{
\t\t\t\t\tnames[k,1:4] = unlist(strsplit(sequencesID[k], "_"))
\t\t\t\t\tsequence = ""
\t\t\t\t\tfor (l in 1:dim(fastaChar)[2]) sequence = paste(sequence,fastaChar[k,l],sep="")
\t\t\t\t\tsequence = gsub("a","A",sequence); sequence = gsub("c","C",sequence)
\t\t\t\t\tsequence = gsub("g","G",sequence); sequence = gsub("t","T",sequence)
\t\t\t\t\tsequences[k,1] = sequence
\t\t\t\t}
\t\t\tIDs = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
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
\t\t\tIDs = names[,1]; years = names[,2]; locations1 = names[,4]; locations2 = names[,3]
\t\t\tcoordinates = cbind(IDs, as.numeric(years), as.numeric(locations1), as.numeric(locations2))
\t\t}
\trasterNames = c()
\tfor (i in 1:length(envVariables))
\t\t{
\t\t\trasterNames = c(rasterNames, names(envVariables[[i]]))
\t\t}
\tfor (i in 1:length(envVariables))
\t\t{
\t\t\tenvVariable = rasterToPolygons(envVariables[[i]])
\t\t\tnames(envVariable) = "envVariable"
\t\t\tif (resistances[[i]] == TRUE) resistance = "R"
\t\t\tif (resistances[[i]] == FALSE) resistance = "C"
\t\t\tcartogram = quick.carto(envVariable, envVariable@data$envVariable, blur=0)
\t\t\tcoords1 = matrix(nrow=dim(coordinates)[1], ncol=2)
\t\t\tfor (j in 1:dim(coordinates)[1])
\t\t\t\t{
\t\t\t\t\tcoords1[j,1] = as.numeric(coordinates[j,3])
\t\t\t\t\tcoords1[j,2] = as.numeric(coordinates[j,4])
\t\t\t\t}
\t\t\tif (showingPlots == TRUE)
\t\t\t\t{
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)[c(1:80)]
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"BuPu"))(100)[c(1:80)]
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"YlGn"))(100)[c(1:80)]
\t\t\t\t\tpts = list("sp.points", SpatialPoints(coords1), pch=19, cex=0.6, col="black")
\t\t\t\t\tspplot(envVariable, border=NA, col.regions=cols, col=NA, sp.layout=list(pts))
\t\t\t\t}
\t\t\tindices = matrix(nrow=dim(coordinates)[1], ncol=1)
\t\t\tfor (j in 1:dim(coords1)[1])
\t\t\t\t{
\t\t\t\t\t# cS = c()
\t\t\t\t\tfound = FALSE
\t\t\t\t\tc1 = cbind(coords1[j,1],coords1[j,2])
\t\t\t\t\twhile (found == FALSE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfor (k in 1:length(envVariable@polygons))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tc2 = envVariable@polygons[[k]]@Polygons[[1]]@coords
\t\t\t\t\t\t\t\t\t# cS = rbind(cS,c)
\t\t\t\t\t\t\t\t\tif (point.in.polygon(c1[1,1],c1[1,2],c2[,1],c2[,2]) == 1)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tindices[j,1] = k; found = TRUE # print(k)
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tif (found == FALSE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tdX = runif(1,-0.001,0.001); dY = runif(1,-0.001,0.001)
\t\t\t\t\t\t\t\t\tc1[1,1] = c1[1,1]+dX; c1[1,2] = c1[1,2]+dY
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t}
\t\t\tcoords2 = matrix(nrow=dim(coordinates)[1], ncol=2)
\t\t\tfor (j in 1:dim(coords1)[1])
\t\t\t\t{
\t\t\t\t\tcoords2[j,] = cartogram@polygons[[indices[j,1]]]@Polygons[[1]]@labpt
\t\t\t\t}
\t\t\tif (showingPlots == TRUE)
\t\t\t\t{
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)[c(1:80)]
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"BuPu"))(100)[c(1:80)]
\t\t\t\t\tcols = colorRampPalette(brewer.pal(9,"YlGn"))(100)[c(1:80)]
\t\t\t\t\tpts = list("sp.points", SpatialPoints(coords2), pch=19, cex=0.6, col="black")
\t\t\t\t\tif (resistances[[i]] == FALSE) cols = rev(cols)
\t\t\t\t\tspplot(cartogram, border=NA, col.regions=cols, col=NA, sp.layout=list(pts))
\t\t\t\t}
\t\t\tcoords3 = coords2
\t\t\tcoords3[,1] = (coords3[,1]-min(coords3[,1]))/(max(coords3[,1])-min(coords3[,1]))
\t\t\tcoords3[,2] = (coords3[,2]-min(coords3[,2]))/(max(coords3[,2])-min(coords3[,2]))
\t\t\tcoords3[,1] = (coords3[,1]*(max(coords1[,1])-min(coords1[,1]))+min(coords1[,1]))
\t\t\tcoords3[,2] = (coords3[,2]*(max(coords1[,2])-min(coords1[,2]))+min(coords1[,2]))
\t\t\tif (textFile == TRUE)
\t\t\t\t{
\t\t\t\t\tsink(file=paste(paste(outputName,"_cartogram_",rasterNames[i],"_",resistance,sep=""),".txt",sep=""))
\t\t\t\t\tfor (j in 1:dim(input)[1])
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste(coords3[j,2],coords3[j,1],sep="\t")); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t\tif (fastaAlignment == TRUE)
\t\t\t\t{
\t\t\t\t\tfasta = input0
\t\t\t\t\tif (length(row.names(fasta)) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfor (j in 1:length(row.names(fasta)))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\trow.names(fasta)[j] = paste(IDs[j],years[j],coords3[j,2],coords3[j,1],sep="_")
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}
\t\t\t\t\tif (length(names(fasta)) > 0)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tfor (j in 1:length(names(fasta)))
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tnames(fasta)[j] = paste(IDs[j],years[j],coords3[j,2],coords3[j,1],sep="_")
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t}\t
\t\t\t\t\tsink(file=paste(paste(outputName,"_cartogram_",rasterNames[i],"_",resistance,sep=""),".fasta",sep=""))
\t\t\t\t\tfor (j in 1:length(row.names(fasta)))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(paste(">",row.names(fasta)[j],sep="")); cat("\\n")
\t\t\t\t\t\t\tcat(sequences[j]); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t\tif (treeFile == TRUE)
\t\t\t\t{
\t\t\t\t\tnew.labels = c()
\t\t\t\t\tfor (j in 1:length(IDs))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tnew.labels = c(new.labels, paste(IDs[j],years[j],coords3[j,2],coords3[j,1],sep="_"))
\t\t\t\t\t\t}
\t\t\t\t\tnewTree = tree
\t\t\t\t\tnewTree$tip.label = new.labels
\t\t\t\t\twrite.tree(newTree, paste(outputName,"_cartogram_",rasterNames[i],"_",resistance,".tree",sep=""))
\t\t\t\t\tnewickTree = scan(paste(outputName,"_cartogram_",rasterNames[i],"_",resistance,".tree",sep=""), what="", sep="\\n", quiet=TRUE)
\t\t\t\t\txmlTemplate = scan("Xml_template.xml", what="", sep="\\n", quiet=TRUE)
\t\t\t\t\txmlBlock = c()
\t\t\t\t\tfor (j in 1:length(IDs))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tsequenceName = paste(IDs[j],years[j],coords3[j,2],coords3[j,1],sep="_")
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t<taxon id=\\"",sequenceName,"\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<date value=\\"",years[j],"\\" direction=\\"forwards\\" units=\\"years\\" precision=\\"0.1\\"/>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lat\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coords3[j,2],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lon\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coords3[j,1],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- START Multivariate diffusion model -->",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"location\\">",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coords3[j,2]," ",coords3[j,1],sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- END Multivariate diffusion model -->",sep=""))
\t\t\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t</taxon>",sep=""))
\t\t\t\t\t\t}
\t\t\t\t\tsink(file=paste(outputName,"_cartogram_",rasterNames[i],"_",resistance,".xml",sep=""))
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
\t\t\t\t\tfor (j in 1:length(xmlTemplate))
\t\t\t\t\t\t{
\t\t\t\t\t\t\tcat(xmlTemplate[j]); cat("\\n")
\t\t\t\t\t\t}
\t\t\t\t\tcat("\\n")
\t\t\t\t\tcat("<beast>"); cat("\\n")
\t\t\t\t\tsink(NULL)
\t\t\t\t}
\t\t}
\tif (treeFile == TRUE)
\t\t{
\t\t\twrite.tree(tree, paste(outputName,"_original.tree",sep=""))
\t\t\tnewickTree = scan(paste(outputName,"_original.tree",sep=""), what="", sep="\\n", quiet=TRUE)
\t\t\txmlTemplate = scan("Xml_template.xml", what="", sep="\\n", quiet=TRUE)
\t\t\txmlBlock = c()
\t\t\tfor (j in 1:length(IDs))
\t\t\t\t{
\t\t\t\t\tsequenceName = paste(coordinates[j,1],coordinates[j,2],coordinates[j,4],coordinates[j,3],sep="_")
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t<taxon id=\\"",sequenceName,"\\">",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<date value=\\"",coordinates[j,2],"\\" direction=\\"forwards\\" units=\\"years\\" precision=\\"0.1\\"/>",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lat\\">",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coordinates[j,4],sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"lon\\">",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coordinates[j,3],sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- START Multivariate diffusion model -->",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<attr name=\\"location\\">",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t\\t",coordinates[j,4]," ",coordinates[j,3],sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t</attr>",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t\\t<!-- END Multivariate diffusion model -->",sep=""))
\t\t\t\t\txmlBlock = c(xmlBlock, paste("\\t</taxon>",sep=""))
\t\t\t\t}
\t\t\tsink(file=paste(outputName,"_noTransformation.xml",sep=""))
\t\t\tcat("<beast>"); cat("\\n")
\t\t\tcat("\\t<taxa id=\\"taxa\\">"); cat("\\n")
\t\t\tfor (j in 1:length(xmlBlock))
\t\t\t\t{
\t\t\t\t\tcat(xmlBlock[j]); cat("\\n")
\t\t\t\t}
\t\t\tcat("\\t</taxa>"); cat("\\n")
\t\t\tcat("\\n")
\t\t\tcat("\\t<newick id=\\"startingTree.sim\\">"); cat("\\n")
\t\t\tcat(newickTree); cat("\\n")
\t\t\tcat("\\t</newick>"); cat("\\n")
\t\t\tcat("\\n")
\t\t\tfor (j in 1:length(xmlTemplate))
\t\t\t\t{
\t\t\t\t\tcat(xmlTemplate[j]); cat("\\n")
\t\t\t\t}
\t\t\tcat("\\n")
\t\t\tcat("<beast>"); cat("\\n")
\t\t\tsink(NULL)
\t\t}\t
}
