cartogramTransformation <-
function (input, envVariables = list(), resistances = c(T), outputName = "") 
{
    input0 = input
    showingPlots = FALSE
    textFile = FALSE
    fastaAlignment = FALSE
    treeFile = FALSE
    for (i in 1:length(envVariables)) {
        if (resistances[[i]] == FALSE) {
            envVariables[[i]][!is.na(envVariables[[i]][])] = 1/envVariables[[i]][!is.na(envVariables[[i]][])]
        }
    }
    if (attr(input, "class") == "data.frame") 
        textFile = TRUE
    if (attr(input, "class") == "DNAbin") 
        fastaAlignment = TRUE
    if (attr(input, "class") == "phylo") 
        treeFile = TRUE
    if (textFile == TRUE) {
        sequencesID = seq(1, dim(input)[1], 1)
        years = rep(0, dim(input)[1])
        locations1 = input[, 2]
        locations2 = input[, 1]
        coordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), 
            as.numeric(locations2))
    }
    if (fastaAlignment == TRUE) {
        fastaChar = as.character(input)
        sequencesID = rownames(fastaChar)
        names = matrix(nrow = dim(fastaChar)[1], ncol = 4)
        sequences = matrix(nrow = dim(fastaChar)[1], ncol = 1)
        for (k in 1:dim(fastaChar)[1]) {
            names[k, 1:4] = unlist(strsplit(sequencesID[k], "_"))
            sequence = ""
            for (l in 1:dim(fastaChar)[2]) sequence = paste(sequence, 
                fastaChar[k, l], sep = "")
            sequence = gsub("a", "A", sequence)
            sequence = gsub("c", "C", sequence)
            sequence = gsub("g", "G", sequence)
            sequence = gsub("t", "T", sequence)
            sequences[k, 1] = sequence
        }
        IDs = names[, 1]
        years = names[, 2]
        locations1 = names[, 4]
        locations2 = names[, 3]
        coordinates = cbind(sequencesID, as.numeric(years), as.numeric(locations1), 
            as.numeric(locations2))
    }
    if (treeFile == TRUE) {
        tree = input
        labels = tree$tip.label
        labels = gsub("'", "", labels)
        names = matrix(nrow = length(labels), ncol = 4)
        for (k in 1:length(labels)) {
            names[k, 1:4] = unlist(strsplit(labels[k], "_"))
        }
        IDs = names[, 1]
        years = names[, 2]
        locations1 = names[, 4]
        locations2 = names[, 3]
        coordinates = cbind(IDs, as.numeric(years), as.numeric(locations1), 
            as.numeric(locations2))
    }
    rasterNames = c()
    for (i in 1:length(envVariables)) {
        rasterNames = c(rasterNames, names(envVariables[[i]]))
    }
    for (i in 1:length(envVariables)) {
        envVariable = rasterToPolygons(envVariables[[i]])
        names(envVariable) = "envVariable"
        if (resistances[[i]] == TRUE) 
            resistance = "R"
        if (resistances[[i]] == FALSE) 
            resistance = "C"
        cartogram = quick.carto(envVariable, envVariable@data$envVariable, 
            blur = 0)
        coords1 = matrix(nrow = dim(coordinates)[1], ncol = 2)
        for (j in 1:dim(coordinates)[1]) {
            coords1[j, 1] = as.numeric(coordinates[j, 3])
            coords1[j, 2] = as.numeric(coordinates[j, 4])
        }
        if (showingPlots == TRUE) {
            cols = colorRampPalette(brewer.pal(9, "YlOrBr"))(100)[c(1:80)]
            cols = colorRampPalette(brewer.pal(9, "BuPu"))(100)[c(1:80)]
            cols = colorRampPalette(brewer.pal(9, "YlGn"))(100)[c(1:80)]
            pts = list("sp.points", SpatialPoints(coords1), pch = 19, 
                cex = 0.6, col = "black")
            spplot(envVariable, border = NA, col.regions = cols, 
                col = NA, sp.layout = list(pts))
        }
        indices = matrix(nrow = dim(coordinates)[1], ncol = 1)
        for (j in 1:dim(coords1)[1]) {
            found = FALSE
            c1 = cbind(coords1[j, 1], coords1[j, 2])
            while (found == FALSE) {
                for (k in 1:length(envVariable@polygons)) {
                  c2 = envVariable@polygons[[k]]@Polygons[[1]]@coords
                  if (point.in.polygon(c1[1, 1], c1[1, 2], c2[, 
                    1], c2[, 2]) == 1) {
                    indices[j, 1] = k
                    found = TRUE
                  }
                }
                if (found == FALSE) {
                  dX = runif(1, -0.001, 0.001)
                  dY = runif(1, -0.001, 0.001)
                  c1[1, 1] = c1[1, 1] + dX
                  c1[1, 2] = c1[1, 2] + dY
                }
            }
        }
        coords2 = matrix(nrow = dim(coordinates)[1], ncol = 2)
        for (j in 1:dim(coords1)[1]) {
            coords2[j, ] = cartogram@polygons[[indices[j, 1]]]@Polygons[[1]]@labpt
        }
        if (showingPlots == TRUE) {
            cols = colorRampPalette(brewer.pal(9, "YlOrBr"))(100)[c(1:80)]
            cols = colorRampPalette(brewer.pal(9, "BuPu"))(100)[c(1:80)]
            cols = colorRampPalette(brewer.pal(9, "YlGn"))(100)[c(1:80)]
            pts = list("sp.points", SpatialPoints(coords2), pch = 19, 
                cex = 0.6, col = "black")
            if (resistances[[i]] == FALSE) 
                cols = rev(cols)
            spplot(cartogram, border = NA, col.regions = cols, 
                col = NA, sp.layout = list(pts))
        }
        coords3 = coords2
        coords3[, 1] = (coords3[, 1] - min(coords3[, 1]))/(max(coords3[, 
            1]) - min(coords3[, 1]))
        coords3[, 2] = (coords3[, 2] - min(coords3[, 2]))/(max(coords3[, 
            2]) - min(coords3[, 2]))
        coords3[, 1] = (coords3[, 1] * (max(coords1[, 1]) - min(coords1[, 
            1])) + min(coords1[, 1]))
        coords3[, 2] = (coords3[, 2] * (max(coords1[, 2]) - min(coords1[, 
            2])) + min(coords1[, 2]))
        if (textFile == TRUE) {
            sink(file = paste(paste(outputName, "_cartogram_", 
                rasterNames[i], "_", resistance, sep = ""), ".txt", 
                sep = ""))
            for (j in 1:dim(input)[1]) {
                cat(paste(coords3[j, 2], coords3[j, 1], sep = "\t"))
                cat("\n")
            }
            sink(NULL)
        }
        if (fastaAlignment == TRUE) {
            fasta = input0
            if (length(row.names(fasta)) > 0) {
                for (j in 1:length(row.names(fasta))) {
                  row.names(fasta)[j] = paste(IDs[j], years[j], 
                    coords3[j, 2], coords3[j, 1], sep = "_")
                }
            }
            if (length(names(fasta)) > 0) {
                for (j in 1:length(names(fasta))) {
                  names(fasta)[j] = paste(IDs[j], years[j], coords3[j, 
                    2], coords3[j, 1], sep = "_")
                }
            }
            sink(file = paste(paste(outputName, "_cartogram_", 
                rasterNames[i], "_", resistance, sep = ""), ".fasta", 
                sep = ""))
            for (j in 1:length(row.names(fasta))) {
                cat(paste(">", row.names(fasta)[j], sep = ""))
                cat("\n")
                cat(sequences[j])
                cat("\n")
            }
            sink(NULL)
        }
        if (treeFile == TRUE) {
            new.labels = c()
            for (j in 1:length(IDs)) {
                new.labels = c(new.labels, paste(IDs[j], years[j], 
                  coords3[j, 2], coords3[j, 1], sep = "_"))
            }
            newTree = tree
            newTree$tip.label = new.labels
            write.tree(newTree, paste(outputName, "_cartogram_", 
                rasterNames[i], "_", resistance, ".tree", sep = ""))
            newickTree = scan(paste(outputName, "_cartogram_", 
                rasterNames[i], "_", resistance, ".tree", sep = ""), 
                what = "", sep = "\n", quiet = TRUE)
            xmlTemplate = scan("Xml_template.xml", what = "", 
                sep = "\n", quiet = TRUE)
            xmlBlock = c()
            for (j in 1:length(IDs)) {
                sequenceName = paste(IDs[j], years[j], coords3[j, 
                  2], coords3[j, 1], sep = "_")
                xmlBlock = c(xmlBlock, paste("\t<taxon id=\"", 
                  sequenceName, "\">", sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<date value=\"", 
                  years[j], "\" direction=\"forwards\" units=\"years\" precision=\"0.1\"/>", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"lat\">", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t\t", coords3[j, 
                  2], sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"lon\">", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t\t", coords3[j, 
                  1], sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<!-- START Multivariate diffusion model -->", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"location\">", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t\t", coords3[j, 
                  2], " ", coords3[j, 1], sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
                xmlBlock = c(xmlBlock, paste("\t\t<!-- END Multivariate diffusion model -->", 
                  sep = ""))
                xmlBlock = c(xmlBlock, paste("\t</taxon>", sep = ""))
            }
            sink(file = paste(outputName, "_cartogram_", rasterNames[i], 
                "_", resistance, ".xml", sep = ""))
            cat("<beast>")
            cat("\n")
            cat("\t<taxa id=\"taxa\">")
            cat("\n")
            for (j in 1:length(xmlBlock)) {
                cat(xmlBlock[j])
                cat("\n")
            }
            cat("\t</taxa>")
            cat("\n")
            cat("\n")
            cat("\t<newick id=\"startingTree.sim\">")
            cat("\n")
            cat(newickTree)
            cat("\n")
            cat("\t</newick>")
            cat("\n")
            cat("\n")
            for (j in 1:length(xmlTemplate)) {
                cat(xmlTemplate[j])
                cat("\n")
            }
            cat("\n")
            cat("<beast>")
            cat("\n")
            sink(NULL)
        }
    }
    if (treeFile == TRUE) {
        write.tree(tree, paste(outputName, "_original.tree", 
            sep = ""))
        newickTree = scan(paste(outputName, "_original.tree", 
            sep = ""), what = "", sep = "\n", quiet = TRUE)
        xmlTemplate = scan("Xml_template.xml", what = "", sep = "\n", 
            quiet = TRUE)
        xmlBlock = c()
        for (j in 1:length(IDs)) {
            sequenceName = paste(coordinates[j, 1], coordinates[j, 
                2], coordinates[j, 4], coordinates[j, 3], sep = "_")
            xmlBlock = c(xmlBlock, paste("\t<taxon id=\"", sequenceName, 
                "\">", sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<date value=\"", 
                coordinates[j, 2], "\" direction=\"forwards\" units=\"years\" precision=\"0.1\"/>", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"lat\">", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t\t", coordinates[j, 
                4], sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"lon\">", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t\t", coordinates[j, 
                3], sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<!-- START Multivariate diffusion model -->", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<attr name=\"location\">", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t\t", coordinates[j, 
                4], " ", coordinates[j, 3], sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t</attr>", sep = ""))
            xmlBlock = c(xmlBlock, paste("\t\t<!-- END Multivariate diffusion model -->", 
                sep = ""))
            xmlBlock = c(xmlBlock, paste("\t</taxon>", sep = ""))
        }
        sink(file = paste(outputName, "_noTransformation.xml", 
            sep = ""))
        cat("<beast>")
        cat("\n")
        cat("\t<taxa id=\"taxa\">")
        cat("\n")
        for (j in 1:length(xmlBlock)) {
            cat(xmlBlock[j])
            cat("\n")
        }
        cat("\t</taxa>")
        cat("\n")
        cat("\n")
        cat("\t<newick id=\"startingTree.sim\">")
        cat("\n")
        cat(newickTree)
        cat("\n")
        cat("\t</newick>")
        cat("\n")
        cat("\n")
        for (j in 1:length(xmlTemplate)) {
            cat(xmlTemplate[j])
            cat("\n")
        }
        cat("\n")
        cat("<beast>")
        cat("\n")
        sink(NULL)
    }
}
