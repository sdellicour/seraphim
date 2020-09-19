torusRandomisation <-
function(rast)\t{

\tnberOfRows = dim(rast)[1]; nberOfColumns = dim(rast)[2]
\tjumpOfRows = as.integer(runif (1, min = 1, nberOfRows))
\tjumpOfColumns = as.integer(runif (1, min = 1, nberOfColumns))
\treflectLines = FALSE; random = runif(1)
\tif (random < 0.5)
\t\t{
\t\t\treflectLines = TRUE
\t\t}
\treflectColumns = FALSE; random = runif(1)
\tif (random < 0.5)
\t\t{
\t\t\treflectColumns = TRUE
\t\t}\t\t\t\t\t\t
\tmat = matrix(0, nrow=dim(rast)[1], ncol=dim(rast)[2])
\tmat_transit = matrix(0, nrow=dim(rast)[1], ncol=dim(rast)[2])
\tmat = raster::as.matrix(rast)
\tmat_transit = mat
\trandom = runif(1); randomTot = FALSE
\trandomColumns = FALSE; randomRows = FALSE
\tcounter = 0
\twhile (randomTot == FALSE)
\t\t{\t
\t\t\tif ((random < 0.5) & (randomColumns == FALSE))
\t\t\t\t{
\t\t\t\t\tcounter = counter + 1
\t\t\t\t\tfor (i in 1:nberOfRows)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tNAs = 0
\t\t\t\t\t\t\tfor (j in 1:nberOfColumns)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (!is.na(mat[i,j]))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tnewJ = j + jumpOfColumns + NAs
\t\t\t\t\t\t\t\t\t\t\tif ((j + jumpOfColumns + NAs) > nberOfColumns)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = j + jumpOfColumns + NAs - nberOfColumns
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\twhile (is.na(mat_transit[i,newJ]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs + 1
\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = newJ + 1
\t\t\t\t\t\t\t\t\t\t\t\t\tif (newJ > nberOfColumns)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = newJ - nberOfColumns
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tif (newJ < 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = newJ + nberOfColumns
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tnewJ = j + jumpOfColumns + NAs
\t\t\t\t\t\t\t\t\t\t\tif ((j + jumpOfColumns + NAs) > nberOfColumns)
\t\t\t\t\t\t\t\t\t\t\t\t{ 
\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = j + jumpOfColumns + NAs - nberOfColumns
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tif ((j + jumpOfColumns + NAs) < 0)
\t\t\t\t\t\t\t\t\t\t\t\t{ 
\t\t\t\t\t\t\t\t\t\t\t\t\tnewJ = j + jumpOfColumns + NAs + nberOfColumns
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tmat[i,j] = mat_transit[i,newJ]
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tNAs = NAs - 1
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tif (reflectLines == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tNAs = 0\t
\t\t\t\t\t\t\t\t\tline_transit = matrix(nrow=1, ncol=dim(rast)[2])
\t\t\t\t\t\t\t\t\tfor (j in 1:nberOfColumns)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (!is.na(mat[i,j]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tJ = nberOfColumns - (j-1) - NAs
\t\t\t\t\t\t\t\t\t\t\t\t\twhile(is.na(mat[i,J]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tJ = nberOfColumns - (j-1) - NAs
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tline_transit[1,j] = mat[i,J]
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs - 1
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tmat[i,] = line_transit
\t\t\t\t\t\t\t\t}\t\t
\t\t\t\t\t\t}
\t\t\t\t\tmat_transit = mat
\t\t\t\t\trandomColumns = TRUE
\t\t\t\t\trandom = 1
\t\t\t\t}\telse if (randomRows == FALSE)\t{
\t\t\t\t\tcounter = counter + 1\t\t\t\t\t\t\t\t\t\t\t
\t\t\t\t\tfor (j in 1:nberOfColumns)
\t\t\t\t\t\t{
\t\t\t\t\t\t\tNAs = 0
\t\t\t\t\t\t\tfor (i in 1:nberOfRows)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tif (!is.na(mat[i,j]))
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tnewI = i + jumpOfRows + NAs
\t\t\t\t\t\t\t\t\t\t\tif ((i + jumpOfRows + NAs) > nberOfRows)
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = i + jumpOfRows + NAs - nberOfRows
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\twhile (is.na(mat_transit[newI,j]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs + 1
\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = newI + 1
\t\t\t\t\t\t\t\t\t\t\t\t\tif (newI > nberOfRows)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = newI - nberOfRows
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tif (newI < 0)
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = newI + nberOfRows
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tnewI = i + jumpOfRows + NAs
\t\t\t\t\t\t\t\t\t\t\tif ((i + jumpOfRows + NAs) > nberOfRows)
\t\t\t\t\t\t\t\t\t\t\t\t{ 
\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = i + jumpOfRows + NAs - nberOfRows
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\tif ((i + jumpOfRows + NAs) < 0)
\t\t\t\t\t\t\t\t\t\t\t\t{ 
\t\t\t\t\t\t\t\t\t\t\t\t\tnewI = i + jumpOfRows + NAs + nberOfRows
\t\t\t\t\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t\t\t\t\t\tmat[i,j] = mat_transit[newI,j]
\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\tNAs = NAs - 1
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\tif (reflectColumns == TRUE)
\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\tNAs = 0\t
\t\t\t\t\t\t\t\t\tcolumn_transit = matrix(nrow=dim(rast)[1], ncol=1)
\t\t\t\t\t\t\t\t\tfor (i in 1:nberOfRows)
\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\tif (!is.na(mat[i,j]))
\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tI = nberOfRows - (i-1) - NAs
\t\t\t\t\t\t\t\t\t\t\t\t\twhile(is.na(mat[I,j]))
\t\t\t\t\t\t\t\t\t\t\t\t\t\t{
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs + 1
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tI = nberOfRows - (i-1) - NAs
\t\t\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t\t\t\tcolumn_transit[i,1] = mat[I,j]
\t\t\t\t\t\t\t\t\t\t\t\t}\telse\t{
\t\t\t\t\t\t\t\t\t\t\t\t\tNAs = NAs - 1
\t\t\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\t\t}
\t\t\t\t\t\t\t\t\tmat[,j] = column_transit
\t\t\t\t\t\t\t\t}\t
\t\t\t\t\t\t}
\t\t\t\t\tmat_transit = mat
\t\t\t\t\trandomRows = TRUE
\t\t\t\t\trandom = 0
\t\t\t\t}
\t\t\tif (randomColumns == TRUE)
\t\t\t\t{
\t\t\t\t\tif (randomRows == TRUE)
\t\t\t\t\t\t{
\t\t\t\t\t\t\trandomTot = TRUE
\t\t\t\t\t\t}
\t\t\t\t}
\t\t}
\trast[,] = mat
\treturn(rast)
}
