library(geomorph)
source("plotRefToTarget.r")

# Example data
# 2d
data(plethodon) 
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)
# 3d
data(scallops)
Y3d.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)

# Simple function that allows a user to choose a single point and it gets plotted
# A is 3D array of coords

picknplot.shape <- function(A, animate = FALSE, PCaxis = 1, Nfr = 20, ...){
  plot.args <- list(...)
  PCA <- prcomp(two.d.array(A))
  PC <- PCA$x[,1:2]
  plot(PC, asp = 1, pch = 19)
  if (animate == F){
    continue <- "y"
    while(continue == "y"){
      cat("Pick a point in the shape space", "\n")
      picked.pts <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 1))
      cat("Picked point coordinates are:", "\n")
      cat(picked.pts, "\n") 
      preds <- shape.predictor(A, x = PC, Intercept = FALSE, pred1 = picked.pts) 
      if (dim(A)[2]==2) {
        plot.args$M1 <- cbind(mshape(A), 0)
        plot.args$M2 <- cbind(preds$pred1, 0)
        class(plot.args$M2) <- "predshape.k2"
        view3d(phi = 0, fov = 30, interactive = FALSE) 
        do.call(plotRefToTarget, plot.args)
      }
      if (dim(A)[2]==3){
        plot.args$M1 <- mshape(A)
        plot.args$M2 <- preds$pred1
        class(plot.args$M2) <- "predshape.k3"
        if(plot.args$method == "TPS"){view3d(phi = 0, fov = 30, interactive = FALSE)}
        else view3d(phi = 0, fov = 30, interactive = TRUE) 
        do.call(plotRefToTarget,  args = plot.args)
      }
      ans <- readline("Save deformation grid as png file (y/n)? ")
      if(ans=="y") {
        file.name <- readline("Please provide file name for saving deformation grid (without quotes) ")
        rgl.snapshot(file = file.name)
      }
      if(ans=="n"){
        try(rgl.close(), silent=T)
      }
      continue <- readline("Do you want to pick another point (y/n)? ")
    }  
  }
  if (animate == TRUE){
    pc.gradient <- seq(min(PC[,PCaxis]), max(PC[,PCaxis]), length.out = Nfr)
    preds <- lapply(pc.gradient, function(y){shape.predictor(A, x = PC[,PCaxis], pred = y)})
    if (dim(A)[2]==2) {
      f1 <- function(x) {
        plot.args$M1 <- cbind(mshape(A), 0)
        plot.args$M2 <- cbind(x$pred, 0)
        view3d(phi = 0, fov = 30, interactive = FALSE) 
        do.call(plotRefToTarget, plot.args)
      }
      try(play3d(lapply(preds, f1), duration = Nfr*0.15), silent = T)
    }
    # NOT SURE ABOUT OPTIONS TO IMPLEMENT FOR 3D YET
    # Also, the 3d version is quite slow  
    if (dim(A)[2]==3){
      f1 <- function(x) {
        plot.args$M1 <- mshape(A)
        plot.args$M2 <- x$pred
        view3d(phi = 0, fov = 30, interactive = FALSE) 
        do.call(plotRefToTarget, plot.args)
      }
      try(play3d(lapply(preds, f1), duration = Nfr*0.15), silent = T)
    }
  }
}

# try it out for 2d (method = "TPS" not available yet in rgl machine)
picknplot.shape(Y.gpa$coords, animate = T, mag = 5, method = "points") 
picknplot.shape(Y.gpa$coords, method = "TPS", mag = 1, outline=plethodon$outline) 

# try it out for 3d
picknplot.shape(Y3d.gpa$coords, method = "TPS", mag = 1) 
scallinks <- matrix(c(1,rep(2:16, each=2),1), nrow=16, byrow=TRUE)
picknplot.shape(Y3d.gpa$coords, method = "TPS", mag = 1, links = scallinks) 

