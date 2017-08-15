library(geomorph)

# Example data
data(plethodon) 
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)

data(scallops)
Y3d.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)

# Simple function that allows a user to choose a single point and it gets plotted
# A is 3D array of coords

picknplot.shape <- function(A, ...){
  plot.args <- list(...)
  PCA <- prcomp(two.d.array(A))
  PC <- PCA$x[,1:2]
  plot(PC, asp = 1, pch = 19)
  cat("Pick a point in the shape space")
  picked.pts <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 0.5))
  preds <- shape.predictor(A, x = PC, Intercept = FALSE, pred1 = picked.pts) 
  k <- dim(A)[2]
  if (k==2) {
    plot.args$M1 <- cbind(mshape(A), 0)
    plot.args$M2 <- cbind(preds$pred1, 0)
    view3d(phi = 0, fov = 30, interactive = FALSE) 
    do.call(plotRefToTarget, plot.args)
  }
  if (k==3){
    plot.args$M1 <- mshape(A)
    plot.args$M2 <- preds$pred1
    view3d(phi = 0, fov = 30, interactive = TRUE) 
    do.call(plotRefToTarget,  args = plot.args)
  }
}

# try it out for 2d (method = "TPS" not available yet in rgl machine)
picknplot.shape(Y.gpa$coords, method = "points", mag = 20) 

# try it out for 3d
picknplot.shape(Y3d.gpa$coords, method = "points", mag = 10) 

