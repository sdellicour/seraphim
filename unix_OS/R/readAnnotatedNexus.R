readAnnotatedNexus <-
function(file, tree.names=NULL) {

\tstrip_annotations = function(text) { 
\t\tannotations = list()
\t\tend = 1
\t\t# Merge node and branch annotations:
\t\ttext = gsub("\\\\[&(.*?)\\\\]:\\\\[&(.*?)\\\\]", ":\\\\[&\\\\1,\\\\2\\\\]", text)
\t\ttext = gsub("\\\\[&(.*?)\\\\]:", ":\\\\[&\\\\1]", text)
\t\tpattern = "\\\\[&.*?\\\\]"
\t\trepeat
\t\t\t{
\t\t\t\tmatch = regexpr(pattern=pattern,text=text)
\t\t\t\tif (!(match[1] > 0))
\t\t\t\t\t{
\t\t\t\t\t\tbreak
\t\t\t\t\t}
\t\t\t\tannotations[[end]] = regmatches(text, match)
\t\t\t\ttext = sub(pattern,paste("[",end,"]",sep=""), text)
\t\t\t\tend = end + 1
\t\t\t}
\t\treturn(list(annotations=annotations,tree=text))
\t}
\t
\tsplit_tree_names = function(text) {   
\t\ttext = gsub(pattern="\\\\[.*?\\\\]=", x=text, replacement="")
\t\ttext = gsub(pattern="^tree", x=text, replacement="")
\t\treturn(text)
\t}
\t
\tsplit_tree_traits = function(text) {
\t\t# Pull out annotation:
\t\ttext = regmatches(text,regexpr(pattern="\\\\[.*?\\\\]",text))
\t\t# Remove leading and trailing delimitors:
\t\ttext = substring(text,3,nchar(text)-1)
\t\treturn(text)
\t}
\t
\tparse_value = function(text) {
\t\tvalue = text
\t\tif (length(grep("^\\\\{",value)))
\t\t\t{
\t\t\t\tsave = value
\t\t\t\tvalue = substring(value, 2, nchar(value)-1)
\t\t\t\tdepth = 0
\t\t\t\tr = regexpr(pattern="\\\\{+",value,perl=TRUE)
\t\t\t\tmatch.length = attr(r, "match.length")
\t\t\t\tif (match.length > 0) depth = match.length
\t\t\t\tif (depth == 0)
\t\t\t\t\t{
\t\t\t\t\t\tsplit = ","
\t\t\t\t\t}\telse\t{
\t\t\t    \t\tsplit = paste("(?<=",rep("\\\\}",depth),")",",","(?=" ,rep("\\\\{",depth),")",sep="")
\t\t\t\t\t}
\t\t\t\tif (depth >= 1)
\t\t\t\t\t{
\t\t\t   \t\t\treturn(save)
\t\t\t\t\t}
\t\t\t\tpart = strsplit(value, split, perl=TRUE)[[1]]
\t\t\t\tvalue = list()
\t\t\t\tfor (i in 1:length(part))
\t\t\t\t\t{
\t\t\t    \t\tvalue[[i]] = parse_value(part[i])
\t\t\t\t\t}
\t\t\t}\telse\t{
\t\t\t\tif (!is.na(suppressWarnings(as.numeric(value))))
\t\t\t\t\t{
\t\t\t\t\t\tvalue = as.numeric(value)
\t\t\t\t\t}\telse\t{
\t\t    \t \t\tvalue = gsub("\\\\\\"","", value)
\t\t    \t\t}
\t\t\t}
\t\treturn(value)
\t}
\t
\tparse_traits = function(text, header=F) {
\t\tif (header == TRUE) text = substring(text,3,nchar(text)-1)
\t\tpattern = "(\\"[^\\"]*\\"+|[^,=\\\\s]+)\\\\s*(=\\\\s*(\\\\{[^=]*\\\\}|\\"[^\\"]*\\"+|[^,]+))?"
\t\trgx = gregexpr(pattern,text,perl=TRUE)
\t\ttraits = list()
\t\tn = length(attr(rgx[[1]],"match.length"))
\t\tstart = attr(rgx[[1]],"capture.start")
\t\tnames = attr(rgx[[1]],"capture.names")
\t\tlength = attr(rgx[[1]],"capture.length")
\t\tnames = attr(rgx[[1]],"capture.names")
\t\tfor (i in 1:n)
\t\t\t{
\t\t\t\ts = start[i,3]
\t\t\t\te = s + length[i,3] - 1
\t\t\t\tvalue = substring(text,s,e)
\t\t\t\ts = start[i,1]
\t\t\t\te = s + length[i,1] - 1
\t\t\t\tkey = substring(text,s,e)
\t\t\t\ttraits[[key]] = parse_value(value)
\t\t\t}
\t\treturn(traits)
\t}
\t
\tannotated_clado_build = function(tp) {
\t    stop(paste("Annotated clado.build is not yet implemented.\\n"))
\t}
\t
\tannotated_tree_build = function(tp) {
\t\tadd.internal = function() {
\t\t\t\tedge[j, 1] <<- current.node
\t\t\t\tedge[j, 2] <<- current.node <<- node <<- node + 1L
\t\t\t\tindex[node] <<- j
\t\t\t\tj <<- j + 1L
\t\t\t}
\t\tadd.terminal = function() {
\t\t\t\tedge[j, 1] <<- current.node
\t\t\t\tedge[j, 2] <<- tip
\t\t\t\tindex[tip] <<- j
\t\t\t\tX = unlist(strsplit(new.tpc[k], ":"))
\t\t\t\ttip.label[tip] <<- X[1]
\t\t\t\tindex = length(X)
\t\t\t\tedge.length[j] <<- as.numeric(X[index])
\t\t\t\tif (length(annotations) > 0)
\t\t\t\t\t{
\t\t\t\t\t\tpermute[[j]] <<- annotations[[as.numeric(X[2])]]# permute traits
\t\t\t\t\t}
\t\t\t\tk <<- k + 1L
\t\t\t\ttip <<- tip + 1L
\t\t\t\tj <<- j + 1L
\t\t\t}
\t\tgo.down = function() {
\t\t\t\tl = index[current.node]
\t\t\t\tX = unlist(strsplit(new.tpc[k], ":"))
\t\t\t\tnode.label[current.node - nb.tip] <<- X[1]
\t\t\t\tindex = length(X)
\t\t\t\tedge.length[l] <<- as.numeric(X[index])
\t\t\t\tif (length(annotations) >  0)
\t\t\t\t\t{
\t\t\t\t\t\tpermute[[l]] <<- annotations[[as.numeric(X[2])]]# permute traits
\t\t\t\t\t}
\t\t\t\tk <<- k + 1L
\t\t\t\tcurrent.node <<- edge[l, 1]
\t\t\t}
\t\tif (!length(grep(",", tp))) {
\t\t\t\tobj = list(edge = matrix(c(2L, 1L), 1, 2))
\t\t\t\ttp = unlist(strsplit(tp, "[\\\\(\\\\):;]"))
\t\t\t\tobj$edge.length = as.numeric(tp[3])
\t\t\t\tobj$Nnode = 1L
\t\t\t\tobj$tip.label = tp[2]
\t\t\t\tif (tp[4] != "") obj$node.label = tp[4]
\t\t\t\tclass(obj) = "phylo"
\t\t\t\treturn(obj)
\t\t\t}
\t\tresult = strip_annotations(tp)
\t\tannotations = result$annotations
\t\tnew.tp.stripped = result$tree
\t\t\t# patched for 0.0 root branch length from BEAST2 (not confirmed)
\t\tnew.tp.stripped = gsub("\\\\]0.0;", "\\\\];", new.tp.stripped)
\t\troot.annotation.number = NULL
\t\tm = regexpr("\\\\[\\\\d+\\\\];", new.tp.stripped)
\t\tif (m != -1)
\t\t\t{
\t\t\t\troot.annotation.number = as.numeric(gsub("\\\\[(\\\\d+)\\\\];", "\\\\1", regmatches(new.tp.stripped, m)))
\t\t\t}
\t\tannotations = lapply(annotations, parse_traits, header=T)
\t\ttp.stripped = gsub("\\\\[.*?\\\\]","",tp)
\t\ttpc = unlist(strsplit(tp.stripped, "[\\\\(\\\\),;]"))
\t\ttpc = tpc[nzchar(tpc)]
\t\tnew.tp.stripped = gsub("\\\\[\\\\d+\\\\];", ";", new.tp.stripped)
\t\tnew.tp.stripped = gsub("\\\\[(\\\\d+)\\\\]","\\\\1:", new.tp.stripped)
\t\tnew.tpc = unlist(strsplit(new.tp.stripped, "[\\\\(\\\\),;]"))
\t\tnew.tpc = new.tpc[nzchar(new.tpc)]\t
\t\ttsp = unlist(strsplit(tp.stripped, NULL))
\t\tskeleton = tsp[tsp %in% c("(", ")", ",", ";")]
\t\tnsk = length(skeleton)
\t\tnb.node = sum(skeleton == ")")
\t\tnb.tip = sum(skeleton == ",") + 1
\t\tnb.edge = nb.node + nb.tip\t
\t\tnode.label = character(nb.node)
\t\ttip.label = character(nb.tip)
\t\tedge.length = numeric(nb.edge)
\t\tedge = matrix(0L, nb.edge, 2)
\t\tcurrent.node = node = as.integer(nb.tip + 1)
\t\tedge[nb.edge, 2] = node
\t\tindex = numeric(nb.edge + 1)
\t\tindex[node] = nb.edge
\t\tj = k = tip = 1L
\t\tpermute = list()
\t\tfor (i in 2:nsk)
\t\t\t{
\t\t    \tif (skeleton[i] == "(")
\t\t    \t\t{
\t\t        \t\tadd.internal()
\t\t\t\t    }
\t\t    \tif (skeleton[i] == ",")
\t\t    \t\t{
\t\t        \t\tif (skeleton[i - 1] != ")")
\t\t        \t\t\t{
\t\t\t\t            \tadd.terminal()
\t\t    \t\t\t    }
\t\t    \t\t}
\t\t\t    if (skeleton[i] == ")")
\t\t\t    \t{
\t\t\t\t        if (skeleton[i - 1] == ",")
\t\t\t\t        \t{
\t\t\t\t\t            add.terminal()
\t\t        \t\t\t    go.down()
\t\t\t\t\t        }
\t\t        \t\tif (skeleton[i - 1] == ")")
\t\t        \t\t\t{
\t\t\t\t\t            go.down()
\t\t\t\t\t        }
\t\t    \t\t}
\t\t\t}
\t\tedge = edge[-nb.edge,]
\t\tobj = list(edge = edge, Nnode = nb.node, tip.label = tip.label)
\t\troot.edge = edge.length[nb.edge]
\t\tedge.length = edge.length[-nb.edge]
\t\tif (!all(is.na(edge.length))) obj$edge.length = edge.length
\t\tif (is.na(node.label[1])) node.label[1] = ""
\t\tif (any(nzchar(node.label))) obj$node.label = node.label
\t\tif (!is.na(root.edge)) obj$root.edge = root.edge
\t\tclass(obj) = "phylo"
\t\tattr(obj, "order") = "cladewise"
\t\tif (!is.null(root.annotation.number))
\t\t\t{
\t\t\t\tobj$root.annotation = annotations[[root.annotation.number]]
\t\t\t}
\t\tobj$annotations = permute
\t\tobj
\t}
\t
\treadAnnotatedTree = function(file="", text=NULL, tree.names=NULL, skip=0, comment.char="#", keep.multi=FALSE, ...) {
\t\tunname = function(treetext)
\t\t\t{
\t\t    \tnc = nchar(treetext)
\t\t    \ttstart = 1
\t\t    \twhile (substr(treetext, tstart, tstart) != "(" && tstart <= nc) tstart = tstart + 1
\t\t    \tif (tstart > 1) return(c(substr(treetext, 1, tstart - 1), substr(treetext, tstart, nc)))
\t\t\t\treturn(c("", treetext))
\t\t\t}
\t\tif (!is.null(text))
\t\t\t{
\t\t    \tif (!is.character(text)) stop("argument `text' must be of mode character")
\t\t    \ttree = text
\t\t\t}\telse\t{
\t\t\t\ttree = scan(file=file, what="", sep="\\n", quiet=T, skip=skip, comment.char=comment.char, ...)
\t\t\t}
\t\tif (identical(tree, character(0)))
\t\t\t{
\t\t\t\twarning("empty character string.")
\t\t\t\treturn(NULL)
\t\t\t}
\t\ttree = gsub("[ \\n\\t]", "", tree)
\t\ttree = gsub("\\\\[&R\\\\]", "", tree)
\t\ttree = unlist(strsplit(tree, NULL))
\t\ty = which(tree == ";")
\t\tNtree = length(y)
\t\tx = c(1, y[-Ntree] + 1)
\t\tif (is.na(y[1])) return(NULL)
\t\tSTRING = character(Ntree)
\t\tfor (i in 1:Ntree) STRING[i] = paste(tree[x[i]:y[i]], sep = "", collapse = "")
\t\ttmp = unlist(lapply(STRING, unname))
\t\ttmpnames = tmp[c(TRUE, FALSE)]
\t\tSTRING = tmp[c(FALSE, TRUE)]
\t\tif (is.null(tree.names) && any(nzchar(tmpnames))) tree.names = tmpnames
\t\tcolon = grep(":", STRING)
\t\tif (!is.null(tree.names))
\t\t\t{
\t\t    \ttraits.text = lapply(tree.names, split_tree_traits)
\t\t\t\ttree.names = lapply(tree.names, split_tree_names)
\t\t\t\ttree.traits = lapply(traits.text, parse_traits)
\t\t\t}
\t\tif (!length(colon))
\t\t\t{
\t\t    \tstop(paste("Annotated clado.build is not yet implemented.\\n"))
\t\t\t\tobj = lapply(STRING, annotated_clado_build)
\t\t\t}\telse if (length(colon) == Ntree)\t{
\t\t    \tobj = lapply(STRING, annotated_tree_build)
\t\t\t}\telse\t{
\t\t    \tobj = vector("list", Ntree)
\t\t   \t\tobj[colon] = lapply(STRING[colon], annotated_tree_build)
\t\t\t\tnocolon = (1:Ntree)[!1:Ntree %in% colon]
\t\t\t\tobj[nocolon] = lapply(STRING[nocolon], function(e) ape::read.tree(text = e))
\t\t\t}
\t\tfor (i in 1:Ntree)
\t\t\t{
\t\t    \tROOT = length(obj[[i]]$tip.label) + 1
\t\t\t    if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1)
\t\t        stop(paste("The tree has apparently singleton node(s): cannot read tree file.\\n  Reading Newick file aborted at tree no.",i))
\t\t\t}
\t\tif (Ntree == 1 && !keep.multi)
\t\t\t{
\t\t\t\tobj = obj[[1]]
\t\t\t}\telse\t{
\t\t    \tif (!is.null(tree.names))
\t\t    \t\t{
\t\t        \t\tnames(obj) = tree.names
\t\t\t\t    }
\t\t    \tclass(obj) = "multiPhylo"
\t\t\t}
\t\tobj
\t}

\tX = scan(file=file, what="", sep="\\n", quiet=T)
\tLEFT = grep("\\\\[", X)
\tRIGHT = grep("\\\\]", X)
\tendblock = grep("END;|ENDBLOCK;", X, ignore.case=T)
\tsemico = grep(";", X)
\ti1 = grep("BEGIN TREES;", X, ignore.case=T)
\ti2 = grep("TRANSLATE", X, ignore.case=T)
\ttranslation = if (length(i2) == 1 && i2 > i1)
\t\tTRUE
\telse FALSE
\tif (translation)
\t\t{
\t    \tend = semico[semico > i2][1]
\t\t    x = X[(i2 + 1):end]
\t\t    x = unlist(strsplit(x, "[,; \\t]"))
\t    \tx = x[nzchar(x)]
\t\t    TRANS = matrix(x, ncol = 2, byrow=TRUE)
\t\t    TRANS[, 2] = gsub("['\\"]", "", TRANS[, 2])
\t\t\tn = dim(TRANS)[1]
\t\t}
\tstart = if (translation)
\t\tsemico[semico > i2][1] + 1
\telse semico[semico > i1][1]
\tend = endblock[endblock > i1][1] - 1
\ttree = X[start:end]
\trm(X)
\ttree = tree[tree != ""]
\tsemico = grep(";", tree)
\tNtree = length(semico)
\tif (Ntree == 1 && length(tree) > 1)
\t\t{
\t\t\tSTRING = paste(tree, collapse = "")
\t\t}\telse\t{
\t    \tif (any(diff(semico) != 1))
\t    \t\t{
\t\t\t        STRING = character(Ntree)
\t    \t\t    s = c(1, semico[-Ntree] + 1)
\t\t\t        j = mapply(":", s, semico)
\t    \t\t    if (is.list(j))
\t        \t\t\t{
\t\t\t            \tfor (i in 1:Ntree) STRING[i] = paste(tree[j[[i]]],collapse = "")
\t    \t\t    \t}\telse {
\t\t\t\t\t\t\tfor (i in 1:Ntree) STRING[i] = paste(tree[j[,i]], collapse = "")
\t\t\t        \t}
\t\t\t\t}
\t\t\telse STRING = tree
\t\t}
\trm(tree)
\tSTRING = STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case=T)]
\tNtree = length(STRING)
\tSTRING = gsub("\\\\[&R\\\\]", "", STRING)
\tnms.annontations.trees = sub(" * = *.*", "", STRING)
\tnms.annontations.trees = sub("^ *tree *", "", nms.annontations.trees, ignore.case=T)
\tnms.trees = sub("\\\\s+\\\\[&.*?\\\\]", "", nms.annontations.trees)
\tif (any(nms.trees != nms.annontations.trees))
\t\t{
\t\t\tannotations.trees = sub(".*\\\\[&", "\\\\[&", nms.annontations.trees)
\t\t\tannotations.trees = lapply(annotations.trees, parse_traits, header=TRUE)
\t\t}\telse\t{
\t\t\tannotations.trees = NULL
\t\t}
\tSTRING = sub("^.*? = *", "", STRING)
\tSTRING = gsub("\\\\s", "", STRING)
\tcolon = grep(":", STRING)
\tif (!length(colon))
\t\t{
\t    \tstop("annotated_clado_build is not yet implemented.\\n")
\t\t    trees = lapply(STRING, annotated_clado_build)
\t\t}\telse if (length(colon) == Ntree)\t{
\t\t\ttrees = lapply(STRING, annotated_tree_build)
\t\t}\telse\t{
\t\t\tstop("Unknown error in readAnnotatedNexus.\\n")
\t\t}
\tfor (i in 1:Ntree)
\t\t{
\t\t\ttr = trees[[i]]
\t\t\tif (!translation) n = length(tr$tip.label)
\t\t\tROOT = n + 1
\t    \tif (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 1)
\t    \t\t{
\t\t\t        stop(paste("The tree has apparently singleton node(s): cannot read tree file.\\n  Reading NEXUS file aborted at tree no.",i,sep=""))
\t\t\t\t}
\t\t}
\tif (Ntree == 1)
\t\t{
\t    \ttrees = trees[[1]]
\t\t    if (translation)
\t\t    \t\t{
\t\t        \ttrees$tip.label = TRANS[, 2][as.numeric(trees$tip.label)]
\t\t\t\t}
\t\t}\telse\t{
\t    \tif (!is.null(tree.names)) names(trees) = tree.names
\t\t\tif (translation)
\t\t\t\t{
\t            \tfor (i in 1:Ntree) trees[[i]]$tip.label = TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
\t\t            trees = .compressTipLabel(trees)
\t   \t\t\t}
\t\t\tclass(trees) = "multiPhylo"
\t\t\tif (!all(nms.trees == "")) names(trees) = nms.trees
\t\t}
\t# Add tree-level annotations back on:
\tif (!is.null(annotations.trees))
\t\t{
\t\t\tif (Ntree == 1)
\t\t\t\t{
\t\t\t\t\ttrees$tree.annotations = annotations.trees[[1]]
\t\t\t\t}\telse\t{
\t\t\t\t\tfor (i in 1:Ntree)
\t\t\t\t\t\t{
\t\t\t\t\t\t\ttrees[[i]]$tree.annotations = annotations.trees[[i]]
\t\t\t\t\t\t}
\t\t\t\t}
\t\t}
\ttrees
}
