spreadGL_input = function(tab, mostRecentSamplingDatum) {

	tab1 = tab; tab2 = tab1[,c("length","startYear","endYear","startLat","startLon","endLat","endLon")]
	tab2$id = seq(1,dim(tab2)[1]); tab2$type = "Internal"; tab2[which(!tab1[,"node2"]%in%tab1[,"node1"]),"type"] = "External"
	tab2$startYear = gsub(" UTC","",date_decimal(tab2$startYear)); tab2$endYear = gsub(" UTC","",date_decimal(tab2$endYear))
	tab2 = tab2[,c("id","type","length","startYear","endYear","startLat","startLon","endLat","endLon")]
	colnames(tab2) = c("id","type","length","start_time","end_time","start_lat","start_lon","end_lat","end_lon")
	return(tab2)
}
