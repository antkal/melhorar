# Outline: Make a function that calculates and plots a predicted shape, based on a user click in tangent space.

# Trial comment by Antigoni :)

library(geomorph)
data(plethodon) 
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)

# Example from shape.predictor to get us started
PCA <- plotTangentSpace(Y.gpa$coords)
PC <- PCA$pc.scores[,1:2]
preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
                         pred1 = c(0.045,-0.02), pred2 = c(-0.025,0.06), pred3 = c(-0.06,-0.04)) 
plotRefToTarget(Y.gpa$consensus, preds$pred1)
plotRefToTarget(Y.gpa$consensus, preds$pred2)
plotRefToTarget(Y.gpa$consensus, preds$pred3)

# get PC scores from regular biplot so can use locator()
PCA <- prcomp(two.d.array((Y.gpa$coords)))
PC <- PCA$x[,1:2]
plot(PC, asp=1, pch=19)
picked.pts <- locator(n = 3) # runs, interactive, choose points in plot
picked.pts <- matrix(unlist(picked.pts), ncol=2) # x is first column, y second
# Assuming you chose 3 points:
preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
                         pred1 = picked.pts[1,], pred2 = picked.pts[2,], pred3 = picked.pts[3,]) 
plotRefToTarget(Y.gpa$consensus, preds$pred1)
plotRefToTarget(Y.gpa$consensus, preds$pred2)
plotRefToTarget(Y.gpa$consensus, preds$pred3)

# To figure out: plot shape immediately after first point is picked
# use while()? or stop() when single point is picked, i.e., length(unlist(picked.pts)) ==2


# Simple function that allows a user to choose a single point and it gets plotted
# A is 3D array of coords
picknplot.shape <- function(A){
  PCA <- prcomp(two.d.array(A))
  PC <- PCA$x[,1:2]
  plot(PC, asp=1, pch=19)
  cat("Pick a point in the shape space")
  picked.pts <- unlist(locator(n=1)) 
  preds <- shape.predictor(A, x= PC, Intercept = FALSE, pred1 = picked.pts) 
  plotRefToTarget(mshape(A), preds$pred1)
}
picknplot.shape(Y.gpa$coords) # try it out

# Can we use rgl for all plotting?
# Trying to plot a 2D image in rgl window
view3d(phi = 0, fov = 30, interactive = FALSE) # sets X,Y in view, with Z into screen. And fov removes too much distortion, and interactive stops rotation by mouse click
plot3d(cbind(PC, 0), zlab = "")
# perfect. 

