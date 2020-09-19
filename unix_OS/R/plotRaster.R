plotRaster <-
function(rast, cols=NULL, colNA="gray90", add=FALSE, new=TRUE, addBox=TRUE, addAxes=FALSE, addLegend=FALSE, legendOnly=FALSE) {
\t
\tif (legendOnly == FALSE)
\t\t{ 
\t\t\tif (add == FALSE)
\t\t\t\t{ 
\t\t\t\t\tif (new == TRUE) dev.new(); # dev.new(width=6, height=rast@nrows*(7.5/rast@ncols))
\t\t\t\t\tif (is.null(cols)) cols = rev(colorRampPalette(brewer.pal(11,"RdBu"))(130))[16:115]
\t\t\t\t\tif (addLegend == TRUE) par(mar=c(0,0,0,0), oma=c(1,3.5,2.5,2.5), mgp=c(0,0.4,0), lwd=0.1, bty="o")
\t\t\t\t\tif (addLegend == FALSE) par(mar=c(0,0,0,0), oma=c(3.5,3.5,3.5,0.5), mgp=c(0,0.4,0), lwd=0.1, bty="o")
\t\t\t\t\traster::plot(rast, col=cols, colNA=colNA, box=F, axes=F, legend=F, interpolate=F, useRaster=T)
\t\t\t\t\tif (addAxes == TRUE) 
\t\t\t\t\t\t{
\t\t\t\t\t\t\taxis(1, c(ceiling(xmin(rast)), floor(xmax(rast))), pos=ymin(rast), mgp=c(0,0.4,0), 
\t\t\t\t\t\t\t\t cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
\t\t\t\t\t\t\taxis(2, c(ceiling(ymin(rast)), floor(ymax(rast))), pos=xmin(rast), mgp=c(0,0.7,0), 
\t\t\t\t\t\t\t\t cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
\t\t\t\t\t\t}
\t\t\t\t\tif (addBox == TRUE) rect(xmin(rast), ymin(rast), xmax(rast), ymax(rast), xpd=T, lwd=0.2)
\t\t\t\t}\telse\t{
\t\t\t\t\tplot(rast, col=cols, colNA=colNA, box=F, axes=F, legend=F, interpolate=F, useRaster=T, add=T)
\t\t\t\t}
\t\t}
\tif ((addLegend == TRUE)|(legendOnly == TRUE)) plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.93,0.94,0.3,0.7),
\t\t\t\t \t\t\t\t legend.args=list(text=""), axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.15, tck=-1, col.axis="gray30", line=0, mgp=c(0,0.5,0)))
}
