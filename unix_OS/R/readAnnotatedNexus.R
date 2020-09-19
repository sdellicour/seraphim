readAnnotatedNexus <-
function(file, tree.names=NULL) {

	strip_annotations = function(text) { 
		annotations = list()
		end = 1
		# Merge node and branch annotations:
		text = gsub("\\\\[&(.*?)\\\\]:\\\\[&(.*?)\\\\]", ":\\\\[&\\\\1,\\\\2\\\\]", text)
		text = gsub("\\\\[&(.*?)\\\\]:", ":\\\\[&\\\\1]", text)
		pattern = "\\\\[&.*?\\\\]"
		repeat
			{
				match = regexpr(pattern=pattern,text=text)
				if (!(match[1] > 0))
					{
						break
					}
				annotations[[end]] = regmatches(text, match)
				text = sub(pattern,paste("[",end,"]",sep=""), text)
				end = end + 1
			}
		return(list(annotations=annotations,tree=text))
	}
	
	split_tree_names = function(text) {   
		text = gsub(pattern="\\\\[.*?\\\\]=", x=text, replacement="")
		text = gsub(pattern="^tree", x=text, replacement="")
		return(text)
	}
	
	split_tree_traits = function(text) {
		# Pull out annotation:
		text = regmatches(text,regexpr(pattern="\\\\[.*?\\\\]",text))
		# Remove leading and trailing delimitors:
		text = substring(text,3,nchar(text)-1)
		return(text)
	}
	
	parse_value = function(text) {
		value = text
		if (length(grep("^\\\\{",value)))
			{
				save = value
				value = substring(value, 2, nchar(value)-1)
				depth = 0
				r = regexpr(pattern="\\\\{+",value,perl=TRUE)
				match.length = attr(r, "match.length")
				if (match.length > 0) depth = match.length
				if (depth == 0)
					{
						split = ","
					}	else	{
			    		split = paste("(?<=",rep("\\\\}",depth),")",",","(?=" ,rep("\\\\{",depth),")",sep="")
					}
				if (depth >= 1)
					{
			   			return(save)
					}
				part = strsplit(value, split, perl=TRUE)[[1]]
				value = list()
				for (i in 1:length(part))
					{
			    		value[[i]] = parse_value(part[i])
					}
			}	else	{
				if (!is.na(suppressWarnings(as.numeric(value))))
					{
						value = as.numeric(value)
					}	else	{
		    	 		value = gsub("\\\\\"","", value)
		    		}
			}
		return(value)
	}
	
	parse_traits = function(text, header=F) {
		if (header == TRUE) text = substring(text,3,nchar(text)-1)
		pattern = "(\"[^\"]*\"+|[^,=\\\\s]+)\\\\s*(=\\\\s*(\\\\{[^=]*\\\\}|\"[^\"]*\"+|[^,]+))?"
		rgx = gregexpr(pattern,text,perl=TRUE)
		traits = list()
		n = length(attr(rgx[[1]],"match.length"))
		start = attr(rgx[[1]],"capture.start")
		names = attr(rgx[[1]],"capture.names")
		length = attr(rgx[[1]],"capture.length")
		names = attr(rgx[[1]],"capture.names")
		for (i in 1:n)
			{
				s = start[i,3]
				e = s + length[i,3] - 1
				value = substring(text,s,e)
				s = start[i,1]
				e = s + length[i,1] - 1
				key = substring(text,s,e)
				traits[[key]] = parse_value(value)
			}
		return(traits)
	}
	
	annotated_clado_build = function(tp) {
	    stop(paste("Annotated clado.build is not yet implemented.\\n"))
	}
	
	annotated_tree_build = function(tp) {
		add.internal = function() {
				edge[j, 1] <<- current.node
				edge[j, 2] <<- current.node <<- node <<- node + 1L
				index[node] <<- j
				j <<- j + 1L
			}
		add.terminal = function() {
				edge[j, 1] <<- current.node
				edge[j, 2] <<- tip
				index[tip] <<- j
				X = unlist(strsplit(new.tpc[k], ":"))
				tip.label[tip] <<- X[1]
				index = length(X)
				edge.length[j] <<- as.numeric(X[index])
				if (length(annotations) > 0)
					{
						permute[[j]] <<- annotations[[as.numeric(X[2])]]# permute traits
					}
				k <<- k + 1L
				tip <<- tip + 1L
				j <<- j + 1L
			}
		go.down = function() {
				l = index[current.node]
				X = unlist(strsplit(new.tpc[k], ":"))
				node.label[current.node - nb.tip] <<- X[1]
				index = length(X)
				edge.length[l] <<- as.numeric(X[index])
				if (length(annotations) >  0)
					{
						permute[[l]] <<- annotations[[as.numeric(X[2])]]# permute traits
					}
				k <<- k + 1L
				current.node <<- edge[l, 1]
			}
		if (!length(grep(",", tp))) {
				obj = list(edge = matrix(c(2L, 1L), 1, 2))
				tp = unlist(strsplit(tp, "[\\\\(\\\\):;]"))
				obj$edge.length = as.numeric(tp[3])
				obj$Nnode = 1L
				obj$tip.label = tp[2]
				if (tp[4] != "") obj$node.label = tp[4]
				class(obj) = "phylo"
				return(obj)
			}
		result = strip_annotations(tp)
		annotations = result$annotations
		new.tp.stripped = result$tree
			# patched for 0.0 root branch length from BEAST2 (not confirmed)
		new.tp.stripped = gsub("\\\\]0.0;", "\\\\];", new.tp.stripped)
		root.annotation.number = NULL
		m = regexpr("\\\\[\\\\d+\\\\];", new.tp.stripped)
		if (m != -1)
			{
				root.annotation.number = as.numeric(gsub("\\\\[(\\\\d+)\\\\];", "\\\\1", regmatches(new.tp.stripped, m)))
			}
		annotations = lapply(annotations, parse_traits, header=T)
		tp.stripped = gsub("\\\\[.*?\\\\]","",tp)
		tpc = unlist(strsplit(tp.stripped, "[\\\\(\\\\),;]"))
		tpc = tpc[nzchar(tpc)]
		new.tp.stripped = gsub("\\\\[\\\\d+\\\\];", ";", new.tp.stripped)
		new.tp.stripped = gsub("\\\\[(\\\\d+)\\\\]","\\\\1:", new.tp.stripped)
		new.tpc = unlist(strsplit(new.tp.stripped, "[\\\\(\\\\),;]"))
		new.tpc = new.tpc[nzchar(new.tpc)]	
		tsp = unlist(strsplit(tp.stripped, NULL))
		skeleton = tsp[tsp %in% c("(", ")", ",", ";")]
		nsk = length(skeleton)
		nb.node = sum(skeleton == ")")
		nb.tip = sum(skeleton == ",") + 1
		nb.edge = nb.node + nb.tip	
		node.label = character(nb.node)
		tip.label = character(nb.tip)
		edge.length = numeric(nb.edge)
		edge = matrix(0L, nb.edge, 2)
		current.node = node = as.integer(nb.tip + 1)
		edge[nb.edge, 2] = node
		index = numeric(nb.edge + 1)
		index[node] = nb.edge
		j = k = tip = 1L
		permute = list()
		for (i in 2:nsk)
			{
		    	if (skeleton[i] == "(")
		    		{
		        		add.internal()
				    }
		    	if (skeleton[i] == ",")
		    		{
		        		if (skeleton[i - 1] != ")")
		        			{
				            	add.terminal()
		    			    }
		    		}
			    if (skeleton[i] == ")")
			    	{
				        if (skeleton[i - 1] == ",")
				        	{
					            add.terminal()
		        			    go.down()
					        }
		        		if (skeleton[i - 1] == ")")
		        			{
					            go.down()
					        }
		    		}
			}
		edge = edge[-nb.edge,]
		obj = list(edge = edge, Nnode = nb.node, tip.label = tip.label)
		root.edge = edge.length[nb.edge]
		edge.length = edge.length[-nb.edge]
		if (!all(is.na(edge.length))) obj$edge.length = edge.length
		if (is.na(node.label[1])) node.label[1] = ""
		if (any(nzchar(node.label))) obj$node.label = node.label
		if (!is.na(root.edge)) obj$root.edge = root.edge
		class(obj) = "phylo"
		attr(obj, "order") = "cladewise"
		if (!is.null(root.annotation.number))
			{
				obj$root.annotation = annotations[[root.annotation.number]]
			}
		obj$annotations = permute
		obj
	}
	
	readAnnotatedTree = function(file="", text=NULL, tree.names=NULL, skip=0, comment.char="#", keep.multi=FALSE, ...) {
		unname = function(treetext)
			{
		    	nc = nchar(treetext)
		    	tstart = 1
		    	while (substr(treetext, tstart, tstart) != "(" && tstart <= nc) tstart = tstart + 1
		    	if (tstart > 1) return(c(substr(treetext, 1, tstart - 1), substr(treetext, tstart, nc)))
				return(c("", treetext))
			}
		if (!is.null(text))
			{
		    	if (!is.character(text)) stop("argument `text' must be of mode character")
		    	tree = text
			}	else	{
				tree = scan(file=file, what="", sep="\\n", quiet=T, skip=skip, comment.char=comment.char, ...)
			}
		if (identical(tree, character(0)))
			{
				warning("empty character string.")
				return(NULL)
			}
		tree = gsub("[ \\n\t]", "", tree)
		tree = gsub("\\\\[&R\\\\]", "", tree)
		tree = unlist(strsplit(tree, NULL))
		y = which(tree == ";")
		Ntree = length(y)
		x = c(1, y[-Ntree] + 1)
		if (is.na(y[1])) return(NULL)
		STRING = character(Ntree)
		for (i in 1:Ntree) STRING[i] = paste(tree[x[i]:y[i]], sep = "", collapse = "")
		tmp = unlist(lapply(STRING, unname))
		tmpnames = tmp[c(TRUE, FALSE)]
		STRING = tmp[c(FALSE, TRUE)]
		if (is.null(tree.names) && any(nzchar(tmpnames))) tree.names = tmpnames
		colon = grep(":", STRING)
		if (!is.null(tree.names))
			{
		    	traits.text = lapply(tree.names, split_tree_traits)
				tree.names = lapply(tree.names, split_tree_names)
				tree.traits = lapply(traits.text, parse_traits)
			}
		if (!length(colon))
			{
		    	stop(paste("Annotated clado.build is not yet implemented.\\n"))
				obj = lapply(STRING, annotated_clado_build)
			}	else if (length(colon) == Ntree)	{
		    	obj = lapply(STRING, annotated_tree_build)
			}	else	{
		    	obj = vector("list", Ntree)
		   		obj[colon] = lapply(STRING[colon], annotated_tree_build)
				nocolon = (1:Ntree)[!1:Ntree %in% colon]
				obj[nocolon] = lapply(STRING[nocolon], function(e) ape::read.tree(text = e))
			}
		for (i in 1:Ntree)
			{
		    	ROOT = length(obj[[i]]$tip.label) + 1
			    if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1)
		        stop(paste("The tree has apparently singleton node(s): cannot read tree file.\\n  Reading Newick file aborted at tree no.",i))
			}
		if (Ntree == 1 && !keep.multi)
			{
				obj = obj[[1]]
			}	else	{
		    	if (!is.null(tree.names))
		    		{
		        		names(obj) = tree.names
				    }
		    	class(obj) = "multiPhylo"
			}
		obj
	}

	X = scan(file=file, what="", sep="\\n", quiet=T)
	LEFT = grep("\\\\[", X)
	RIGHT = grep("\\\\]", X)
	endblock = grep("END;|ENDBLOCK;", X, ignore.case=T)
	semico = grep(";", X)
	i1 = grep("BEGIN TREES;", X, ignore.case=T)
	i2 = grep("TRANSLATE", X, ignore.case=T)
	translation = if (length(i2) == 1 && i2 > i1)
		TRUE
	else FALSE
	if (translation)
		{
	    	end = semico[semico > i2][1]
		    x = X[(i2 + 1):end]
		    x = unlist(strsplit(x, "[,; \t]"))
	    	x = x[nzchar(x)]
		    TRANS = matrix(x, ncol = 2, byrow=TRUE)
		    TRANS[, 2] = gsub("['\"]", "", TRANS[, 2])
			n = dim(TRANS)[1]
		}
	start = if (translation)
		semico[semico > i2][1] + 1
	else semico[semico > i1][1]
	end = endblock[endblock > i1][1] - 1
	tree = X[start:end]
	rm(X)
	tree = tree[tree != ""]
	semico = grep(";", tree)
	Ntree = length(semico)
	if (Ntree == 1 && length(tree) > 1)
		{
			STRING = paste(tree, collapse = "")
		}	else	{
	    	if (any(diff(semico) != 1))
	    		{
			        STRING = character(Ntree)
	    		    s = c(1, semico[-Ntree] + 1)
			        j = mapply(":", s, semico)
	    		    if (is.list(j))
	        			{
			            	for (i in 1:Ntree) STRING[i] = paste(tree[j[[i]]],collapse = "")
	    		    	}	else {
							for (i in 1:Ntree) STRING[i] = paste(tree[j[,i]], collapse = "")
			        	}
				}
			else STRING = tree
		}
	rm(tree)
	STRING = STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case=T)]
	Ntree = length(STRING)
	STRING = gsub("\\\\[&R\\\\]", "", STRING)
	nms.annontations.trees = sub(" * = *.*", "", STRING)
	nms.annontations.trees = sub("^ *tree *", "", nms.annontations.trees, ignore.case=T)
	nms.trees = sub("\\\\s+\\\\[&.*?\\\\]", "", nms.annontations.trees)
	if (any(nms.trees != nms.annontations.trees))
		{
			annotations.trees = sub(".*\\\\[&", "\\\\[&", nms.annontations.trees)
			annotations.trees = lapply(annotations.trees, parse_traits, header=TRUE)
		}	else	{
			annotations.trees = NULL
		}
	STRING = sub("^.*? = *", "", STRING)
	STRING = gsub("\\\\s", "", STRING)
	colon = grep(":", STRING)
	if (!length(colon))
		{
	    	stop("annotated_clado_build is not yet implemented.\\n")
		    trees = lapply(STRING, annotated_clado_build)
		}	else if (length(colon) == Ntree)	{
			trees = lapply(STRING, annotated_tree_build)
		}	else	{
			stop("Unknown error in readAnnotatedNexus.\\n")
		}
	for (i in 1:Ntree)
		{
			tr = trees[[i]]
			if (!translation) n = length(tr$tip.label)
			ROOT = n + 1
	    	if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 1)
	    		{
			        stop(paste("The tree has apparently singleton node(s): cannot read tree file.\\n  Reading NEXUS file aborted at tree no.",i,sep=""))
				}
		}
	if (Ntree == 1)
		{
	    	trees = trees[[1]]
		    if (translation)
		    		{
		        	trees$tip.label = TRANS[, 2][as.numeric(trees$tip.label)]
				}
		}	else	{
	    	if (!is.null(tree.names)) names(trees) = tree.names
			if (translation)
				{
	            	for (i in 1:Ntree) trees[[i]]$tip.label = TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
		            trees = .compressTipLabel(trees)
	   			}
			class(trees) = "multiPhylo"
			if (!all(nms.trees == "")) names(trees) = nms.trees
		}
	# Add tree-level annotations back on:
	if (!is.null(annotations.trees))
		{
			if (Ntree == 1)
				{
					trees$tree.annotations = annotations.trees[[1]]
				}	else	{
					for (i in 1:Ntree)
						{
							trees[[i]]$tree.annotations = annotations.trees[[i]]
						}
				}
		}
	trees
}
