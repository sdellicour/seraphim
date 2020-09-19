vectorRandomisation <-
function(rast, fromCoor, toCoor, showingPlots=FALSE)\t{

\t\tfromCoorRand = fromCoor; toCoorRand = toCoor
\t\tfixedStartLocation = TRUE
\t\trotation = function(pt1, pt2, angle)
\t\t\t{
\t\t\t\ts = sin(angle); c = cos(angle)
\t\t\t\tx = pt2[1]-pt1[1]; y = pt2[2]-pt1[2]
\t\t\t\tx_new = (x*c)-(y*s); y_new = (x*s)+(y*c)
\t\t\t\tx_new = x_new+pt1[1]; y_new = y_new+pt1[2]
\t\t\t\treturn(c(x_new,y_new))
\t\t\t}
\t\tif (showingPlots == TRUE) plotRaster(rast)
\t\tsdX = sd(rbind(fromCoor[,1],toCoor[,1])); sdY = sd(rbind(fromCoor[,2],toCoor[,2]))
\t\trangeX = max(rbind(fromCoor[,1],toCoor[,1]))-min(rbind(fromCoor[,1],toCoor[,1]))
\t\trangeY = max(rbind(fromCoor[,2],toCoor[,2]))-min(rbind(fromCoor[,2],toCoor[,2]))
\t\ttwoPointsOnTheGrid = FALSE
\t\tfromCoorRand[,] = NA; toCoorRand[,] = NA
\t\tfor (i in 1:dim(fromCoorRand)[1])
\t\t\t{
\t\t\t\ttwoPointsOnTheGrid = FALSE
\t\t\t\twhile(twoPointsOnTheGrid == FALSE)
\t\t\t\t\t{
\t\t\t\t\t\tpt01 = c(fromCoor[i,1], fromCoor[i,2])
\t\t\t\t\t\tpt02 = c(toCoor[i,1], toCoor[i,2])
\t\t\t\t\t\txTranslation = pt02[1]-pt01[1]
\t\t\t\t\t\tyTranslation = pt02[2]-pt01[2]
\t\t\t\t\t\tpt1 = c(fromCoor[i,1],fromCoor[i,2])
\t\t\t\t\t\tif (fixedStartLocation == FALSE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tpt1NAarea = TRUE
\t\t\t\t\t\t\t\twhile(pt1NAarea == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tpt1_translated = pt1
\t\t\t\t\t\t\t\t\t\txRand = rnorm(1,0,sdX); yRand = rnorm(1,0,sdY); # print(c(xRand,yRand))
\t\t\t\t\t\t\t\t\t\t# xRand = (runif(1)*rangeX*2)-rangeX; yRand = (runif(1)*rangeY*2)-rangeY
\t\t\t\t\t\t\t\t\t\tpt1_translated[1] = pt1[1]+xRand; pt1_translated[2] = pt1[2]+yRand
\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\tif (!is.na(raster::extract(rast,cbind(pt1_translated[1],pt1_translated[2]))))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tpt1NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t\tpt1 = pt1_translated
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t}
\t\t\t\t\t\tpt2 = c(NA,NA)\t
\t\t\t\t\t\tpt2[1] = pt1[1]+xTranslation
\t\t\t\t\t\tpt2[2] = pt1[2]+yTranslation
\t\t\t\t\t\tcounter = 0
\t\t\t\t\t\tpt2NAarea = TRUE
\t\t\t\t\t\twhile(pt2NAarea == TRUE)
\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\tonTheGrid = TRUE
\t\t\t\t\t\t\t\tangle = (2*pi)*runif(1)
\t\t\t\t\t\t\t\tpt2_rotated = rotation(pt1, pt2, angle)
\t\t\t\t\t\t\t\tif (pt2_rotated[1] > extent(rast)@xmax)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (pt2_rotated[1] < extent(rast)@xmin)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\tif (pt2_rotated[2] > extent(rast)@ymax)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (pt2_rotated[2] < extent(rast)@ymin)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tonTheGrid = FALSE
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (onTheGrid == TRUE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\t\t\t\t\t\tNAarea = FALSE
\t\t\t\t\t\t\t\t\t\tif (is.na(raster::extract(rast,cbind(pt2_rotated[1],pt2_rotated[2]))))
\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\tNAarea = TRUE
\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\tif (NAarea == FALSE)
\t\t\t\t\t\t\t\t\t\t\t{\t
\t\t\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (pt2NAarea == FALSE)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\ttwoPointsOnTheGrid = TRUE
\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\tcounter = counter+1\t\t\t\t\t
\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\tif (counter == 100)
\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\tpt2NAarea = FALSE
\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t}
\t\t\t\t\t}
\t\t\t\tpt2 = pt2_rotated
\t\t\t\tif (showingPlots == TRUE)
\t\t\t\t\t{
\t\t\t\t\t\tpoints(pt1[1], pt1[2], pch=16, col='red', cex=0.6)
\t\t\t\t\t\tpoints(pt2[1], pt2[2], pch=16, col='red', cex=0.6)
\t\t\t\t\t\tsegments(pt1[1], pt1[2], pt2[1], pt2[2], col='red', lwd=1)
\t\t\t\t\t}
\t\t\t\tfromCoorRand[i,1] = pt1[1]; fromCoorRand[i,2] = pt1[2]
\t\t\t\ttoCoorRand[i,1] = pt2[1]; toCoorRand[i,2] = pt2[2]
\t\t\t\td1 = sqrt(((pt01[1]-pt02[1])^2)+((pt01[2]-pt02[2])^2))
\t\t\t\td2 = sqrt(((pt1[1]-pt2[1])^2)+((pt1[2]-pt2[2])^2))
\t\t\t}
\t\tif (showingPlots == TRUE) dev.off()
\t\treturn(cbind(fromCoorRand, toCoorRand))\t\t
\t}
