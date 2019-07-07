vectorRandomisation <-
function (rast, fromCoor, toCoor, showingPlots = FALSE) 
{
    fromCoorRand = fromCoor
    toCoorRand = toCoor
    fixedStartLocation = TRUE
    rotation = function(pt1, pt2, angle) {
        s = sin(angle)
        c = cos(angle)
        x = pt2[1] - pt1[1]
        y = pt2[2] - pt1[2]
        x_new = (x * c) - (y * s)
        y_new = (x * s) + (y * c)
        x_new = x_new + pt1[1]
        y_new = y_new + pt1[2]
        return(c(x_new, y_new))
    }
    if (showingPlots == TRUE) 
        plotRaster(rast)
    sdX = sd(rbind(fromCoor[, 1], toCoor[, 1]))
    sdY = sd(rbind(fromCoor[, 2], toCoor[, 2]))
    rangeX = max(rbind(fromCoor[, 1], toCoor[, 1])) - min(rbind(fromCoor[, 
        1], toCoor[, 1]))
    rangeY = max(rbind(fromCoor[, 2], toCoor[, 2])) - min(rbind(fromCoor[, 
        2], toCoor[, 2]))
    twoPointsOnTheGrid = FALSE
    fromCoorRand[, ] = NA
    toCoorRand[, ] = NA
    for (i in 1:dim(fromCoorRand)[1]) {
        twoPointsOnTheGrid = FALSE
        while (twoPointsOnTheGrid == FALSE) {
            pt01 = c(fromCoor[i, 1], fromCoor[i, 2])
            pt02 = c(toCoor[i, 1], toCoor[i, 2])
            xTranslation = pt02[1] - pt01[1]
            yTranslation = pt02[2] - pt01[2]
            pt1 = c(fromCoor[i, 1], fromCoor[i, 2])
            if (fixedStartLocation == FALSE) {
                pt1NAarea = TRUE
                while (pt1NAarea == TRUE) {
                  pt1_translated = pt1
                  xRand = rnorm(1, 0, sdX)
                  yRand = rnorm(1, 0, sdY)
                  pt1_translated[1] = pt1[1] + xRand
                  pt1_translated[2] = pt1[2] + yRand
                  NAarea = FALSE
                  if (!is.na(extract(rast, cbind(pt1_translated[1], 
                    pt1_translated[2])))) {
                    pt1NAarea = FALSE
                    pt1 = pt1_translated
                  }
                }
            }
            pt2 = c(NA, NA)
            pt2[1] = pt1[1] + xTranslation
            pt2[2] = pt1[2] + yTranslation
            counter = 0
            pt2NAarea = TRUE
            while (pt2NAarea == TRUE) {
                onTheGrid = TRUE
                angle = (2 * pi) * runif(1)
                pt2_rotated = rotation(pt1, pt2, angle)
                if (pt2_rotated[1] > extent(rast)@xmax) {
                  onTheGrid = FALSE
                }
                if (pt2_rotated[1] < extent(rast)@xmin) {
                  onTheGrid = FALSE
                }
                if (pt2_rotated[2] > extent(rast)@ymax) {
                  onTheGrid = FALSE
                }
                if (pt2_rotated[2] < extent(rast)@ymin) {
                  onTheGrid = FALSE
                }
                if (onTheGrid == TRUE) {
                  NAarea = FALSE
                  if (is.na(extract(rast, cbind(pt2_rotated[1], 
                    pt2_rotated[2])))) {
                    NAarea = TRUE
                  }
                  if (NAarea == FALSE) {
                    pt2NAarea = FALSE
                  }
                }
                if (pt2NAarea == FALSE) {
                  twoPointsOnTheGrid = TRUE
                }
                else {
                  counter = counter + 1
                }
                if (counter == 100) {
                  pt2NAarea = FALSE
                }
            }
        }
        pt2 = pt2_rotated
        if (showingPlots == TRUE) {
            points(pt1[1], pt1[2], pch = 16, col = "red", cex = 0.6)
            points(pt2[1], pt2[2], pch = 16, col = "red", cex = 0.6)
            segments(pt1[1], pt1[2], pt2[1], pt2[2], col = "red", 
                lwd = 1)
        }
        fromCoorRand[i, 1] = pt1[1]
        fromCoorRand[i, 2] = pt1[2]
        toCoorRand[i, 1] = pt2[1]
        toCoorRand[i, 2] = pt2[2]
        d1 = sqrt(((pt01[1] - pt02[1])^2) + ((pt01[2] - pt02[2])^2))
        d2 = sqrt(((pt1[1] - pt2[1])^2) + ((pt1[2] - pt2[2])^2))
    }
    if (showingPlots == TRUE) 
        dev.off()
    return(cbind(fromCoorRand, toCoorRand))
}
