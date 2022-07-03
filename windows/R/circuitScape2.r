circuitScape2 = function(envVariable, envVariableName, resistance=TRUE, avgResistance=TRUE, fourCells=FALSE, fromCoor, toCoor, OS="Unix", prefix="", ID="", nberOfCores_CS=1)	{
	
	mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=1)
	folder = paste(prefix, "_CStemp_", ID, sep="")
	if (folder%in%dir(getwd())==FALSE) dir.create(file.path(getwd(), folder))
	cs.text = file(paste(folder, "CS_script.jl", sep="/"))
	line1 = "using Circuitscape"; line2 = paste0("compute(\"",folder,"/CS_temporary.ini\")")
	writeLines(c(line1, line2), cs.text)
	close(cs.text)
	data(cs2) # cs2.ini = scan(file="CS_template2.ini",what="",sep="\n",quiet=T)
	index = which(grepl("point_file = ",cs2.ini))
	cs2.ini[index] = paste("point_file = ", getwd(),paste("/",folder,"/focal_points_temp.txt",sep=""),sep="")
	index = which(grepl("output_file = ",cs2.ini))
	cs2.ini[index] = paste("output_file = ",getwd(),paste("/",folder,"/raster_file_temp.txt",sep=""),sep="")
	index = which(grepl("log_file = XXX",cs2.ini))
	cs2.ini[index] = paste("log_file = ",getwd(),paste("/",folder,"/raster_file_temp.log",sep=""),sep="")
	index = which(grepl("included_pairs_file = XXX",cs2.ini))
	cs2.ini[index] = paste("included_pairs_file = ",getwd(),paste("/",folder,"/pairs_to_include.txt",sep=""),sep="")
	index = which(grepl("max_parallel = ",cs2.ini))
	cs2.ini[index] = paste("max_parallel = ",nberOfCores_CS,sep="")
	index = which(grepl("habitat_file = ",cs2.ini))
	if (grepl(".asc",envVariableName)|grepl(".tif",envVariableName)|grepl(".gri",envVariableName))
		{
			cs2.ini[index] = paste("habitat_file = ",getwd(),"/",envVariableName,sep="")
		}	else	{
			cs2.ini[index] = paste("habitat_file = ",getwd(),"/",envVariableName,".asc",sep="")
		}
	if (OS == "Windows") cs2.ini = gsub("/", "\\", cs2.ini, fixed=T)
	index = which(grepl("habitat_map_is_resistances = ",cs2.ini))
	if (resistance == TRUE)
		{
			cs2.ini[index] = "habitat_map_is_resistances = True" 
		}	else	{
			cs2.ini[index] = "habitat_map_is_resistances = False"
		}
	index = which(grepl("connect_using_avg_resistances = ",cs2.ini))
	if (avgResistance == TRUE)
		{
			cs2.ini[index] = "connect_using_avg_resistances = True"
		}	else	{
			cs2.ini[index] = "connect_using_avg_resistances = False"
		}
	index = which(grepl("connect_four_neighbors_only = ",cs2.ini))
	if (fourCells == TRUE)
		{
			cs2.ini[index] = "connect_four_neighbors_only = True"
		}	else	{
			cs2.ini[index] = "connect_four_neighbors_only = False"
		}
	cs.text = file(paste(folder, "CS_temporary.ini", sep="/"))
	writeLines(cs2.ini, cs.text)
	close(cs.text)
	focal_points = paste(folder, "focal_points_temp.txt", sep="/")
	lines1 = cbind(1:length(fromCoor[,1]), fromCoor)
	lines2 = cbind((length(fromCoor[,1])+1):(length(fromCoor[,1])+length(toCoor[,1])), toCoor)
	lines = rbind(lines1, lines2)
	write.table(lines, focal_points, row.names=F, col.names=F)
	pairs_to_include = paste(folder, "pairs_to_include.txt", sep="/")
	tab = matrix(nrow=dim(mat)[1], ncol=2)
	colnames(tab) = c("mode","include")
	tab[,1] = 1:dim(mat)[1]
	tab[,2] = (dim(mat)[1]+1):(2*dim(mat)[1])
	write.table(tab, pairs_to_include, row.names=F, col.names=T, sep="\t", quote=F)
	if (OS == "Unix")
		{
			system(paste("julia",paste(folder,"CS_script.jl", sep="/"),sep=" "), ignore.stdout=T, ignore.stderr=T)
		}		
	tab = read.table(paste0(folder,"/raster_file_temp.txt_resistances.out"), header=F)
	tab = tab[2:dim(tab)[1], 2:dim(tab)[2]]
	sameCoordinates = FALSE
	if (sum(fromCoor-toCoor) == 0) sameCoordinates = TRUE
	if (sameCoordinates == FALSE)
    	{
			for (i in 1:length(fromCoor[,1]))
				{
					mat[i] = tab[i,(i+length(fromCoor[,1]))]
				}
		}	else	{
			mat = matrix(nrow=dim(as.matrix(fromCoor))[1], ncol=dim(as.matrix(fromCoor))[1])
			mat = tab[1:dim(fromCoor)[1],1:dim(fromCoor)[1]]
		}
	if (OS == "Unix")
		{		
			system(paste0("rm -rf ",folder))
		}
	if (OS == "Windows")
		{
			unlink(folder, recursive=T, force=T)
		}
	return(mat)			
}
