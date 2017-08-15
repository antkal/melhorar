library(geomorph)

# Example data
# 2d
data(plethodon) 
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)
# 3d
data(scallops)
Y3d.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)

# Function that creates an animation runing along the selected PCA axis
# A is 3D array of coords
# One needs to select a PCA axis to visualize (PCaxis)
# Nfr is the number of frames to be sampled across the axis
# Ellipses allow access to plotRefToTarget plotting options

shape.animation <- function(A, PCaxis = 1, Nfr = 20, ...){
  plot.args <- list(...)
  PCA <- prcomp(two.d.array(A))
  PC <- PCA$x
  pc.gradient <- seq(min(PC[,PCaxis]), max(PC[,PCaxis]), length.out = Nfr)
  preds <- lapply(pc.gradient, function(y){shape.predictor(A, x = PC[,PCaxis], pred = y)})
  k <- dim(A)[2]
  if (k==2) {
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
  if (k==3){
  f1 <- function(x) {
    plot.args$M1 <- mshape(A)
    plot.args$M2 <- x$pred
    view3d(phi = 0, fov = 30, interactive = FALSE) 
    do.call(plotRefToTarget, plot.args)
    }
    try(play3d(lapply(preds, f1), duration = Nfr*0.15), silent = T)
  }
}

# Check it out
shape.animation(Y.gpa$coords, PCaxis = 1, Nfr = 20, mag = 5, method = "vector")
shape.animation(Y3d.gpa$coords, PCaxis = 1, Nfr = 20, mag = 5, method = "vector")
