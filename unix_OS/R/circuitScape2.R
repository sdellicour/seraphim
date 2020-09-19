circuitScape2 <-
function(envVariable, envVariableName, resistance=TRUE, avgResistance=TRUE, fourCells=FALSE, fromCoor, toCoor, OS="Unix", prefix="", ID="", nberOfCores_CS=1)\t{
\t
\tmat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=1)
\tfolder = paste(prefix, "_CStemp_", ID, sep="")
\tif (folder%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), folder))
\tcs.text = file(paste(folder, "CS_script.jl", sep="/"))
\tline1 = "using Circuitscape"; line2 = paste0("compute(\\"",folder,"/CS_temporary.ini\\")")
\twriteLines(c(line1, line2), cs.text)
\tclose(cs.text)
\tdata(cs2) # cs2.ini = scan(file="CS_template2.ini",what="",sep="\\n",quiet=T)
\tindex = which(grepl("point_file = ",cs2.ini))
\tcs2.ini[index] = paste("point_file = ", getwd(),paste("/",folder,"/focal_points_temp.txt",sep=""),sep="")
\tindex = which(grepl("output_file = ",cs2.ini))
\tcs2.ini[index] = paste("output_file = ",getwd(),paste("/",folder,"/raster_file_temp.txt",sep=""),sep="")
\tindex = which(grepl("log_file = XXXX",cs2.ini))
\tcs2.ini[index] = paste("log_file = ",getwd(),paste("/",folder,"/raster_file_temp.log",sep=""),sep="")
\tindex = which(grepl("included_pairs_file = XXXX",cs2.ini))
\tcs2.ini[index] = paste("included_pairs_file = ",getwd(),paste("/",folder,"/pairs_to_include.txt",sep=""),sep="")
\tindex = which(grepl("max_parallel = ",cs2.ini))
\tcs2.ini[index] = paste("max_parallel = ",nberOfCores_CS,sep="")
\tindex = which(grepl("habitat_file = ",cs2.ini))
\tif (grepl(".asc",envVariableName)|grepl(".tif",envVariableName)|grepl(".gri",envVariableName))
\t\t{
\t\t\tcs2.ini[index] = paste("habitat_file = ",getwd(),"/",envVariableName,sep="")
\t\t}\telse\t{
\t\t\tcs2.ini[index] = paste("habitat_file = ",getwd(),"/",envVariableName,".asc",sep="")
\t\t}
\tif (OS == "Windows") cs2.ini = gsub("/", "\\\\", cs2.ini, fixed=T)
\tindex = which(grepl("habitat_map_is_resistances = ",cs2.ini))
\tif (resistance == TRUE)
\t\t{
\t\t\t# cs2.ini[index] = "ground_file_is_resistances = False"
\t\t\tcs2.ini[index] = "habitat_map_is_resistances = True" 
\t\t}\telse\t{
\t\t\t# cs2.ini[index] = "ground_file_is_resistances = False"
\t\t\tcs2.ini[index] = "habitat_map_is_resistances = False"
\t\t}
\tindex = which(grepl("connect_using_avg_resistances = ",cs2.ini))
\tif (avgResistance == TRUE)
\t\t{
\t\t\tcs2.ini[index] = "connect_using_avg_resistances = True"
\t\t}\telse\t{
\t\t\tcs2.ini[index] = "connect_using_avg_resistances = False"
\t\t}
\tindex = which(grepl("connect_four_neighbors_only = ",cs2.ini))
\tif (fourCells == TRUE)
\t\t{
\t\t\tcs2.ini[index] = "connect_four_neighbors_only = True"
\t\t}\telse\t{
\t\t\tcs2.ini[index] = "connect_four_neighbors_only = False"
\t\t}
\tcs.text = file(paste(folder, "CS_temporary.ini", sep="/"))
\twriteLines(cs2.ini, cs.text)
\tclose(cs.text)
\tfocal_points = paste(folder, "focal_points_temp.txt", sep="/")
\tlines1 = cbind(1:length(fromCoor[,1]), fromCoor)
\tlines2 = cbind((length(fromCoor[,1])+1):(length(fromCoor[,1])+length(toCoor[,1])), toCoor)
\tlines = rbind(lines1, lines2)
\twrite.table(lines, focal_points, row.names=F, col.names=F)
\tpairs_to_include = paste(folder, "pairs_to_include.txt", sep="/")
\ttab = matrix(nrow=dim(mat)[1], ncol=2)
\tcolnames(tab) = c("mode","include")
\ttab[,1] = 1:dim(mat)[1]
\ttab[,2] = (dim(mat)[1]+1):(2*dim(mat)[1])
\twrite.table(tab, pairs_to_include, row.names=F, col.names=T, sep="\\t", quote=F)
\tif (OS == "Unix")
\t\t{
\t\t\tsystem(paste("julia",paste(folder,"CS_script.jl", sep="/"),sep=" "), ignore.stdout=T, ignore.stderr=T)
\t\t\t# system(paste("python2",paste(folder,"CS_script2.sh", sep="/"),sep=" "), ignore.stdout=T, ignore.stderr=T)
\t\t}\t\t
\ttab = read.table(paste0(folder,"/raster_file_temp_resistances.out"), header=F)
\ttab = tab[2:dim(tab)[1], 2:dim(tab)[2]]
\tsameCoordinates = FALSE
\tif (sum(fromCoor-toCoor) == 0) sameCoordinates = TRUE
\tif (sameCoordinates == FALSE)
    \t{
\t\t\tfor (i in 1:length(fromCoor[,1]))
\t\t\t\t{
\t\t\t\t\tmat[i] = tab[i,(i+length(fromCoor[,1]))]
\t\t\t\t}
\t\t}\telse\t{
\t\t\tmat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=dim(as.matrix(fromCoor))[1])
\t\t\tmat = tab[1:dim(fromCoor)[1],1:dim(fromCoor)[1]]
\t\t}
\tif (OS == "Unix")
\t\t{\t\t
\t\t\tsystem(paste0("rm -rf ",folder))
\t\t}
\tif (OS == "Windows")
\t\t{
\t\t\tunlink(folder, recursive=T, force=T)
\t\t}
\treturn(mat)\t\t\t
}
