# Started setting it up, but I have to run to a meeting and I have not really advanced
# Will come back to it soon ;)

library(geomorph)

# Example data (with phylo)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)


gm.prcomp <- function (A, phy = NULL, ...){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }
  dots <- list(...)
  retx <- dots$retx
  if(is.null(retx)) retx <- TRUE
  scale. <- dots$scale.
  if(is.null(scale.)) scale. <- FALSE
  center <- dots$center
  if(is.null(center)) center <- TRUE
  tol <- dots$tol
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  if(is.null(tol)){
    d <- prcomp(x)$sdev^2
    cd <-cumsum(d)/sum(d)
    cd <- length(which(cd < 1)) 
    if(length(cd) < length(d)) cd <- cd + 1
    tol <- max(c(d[cd]/d[1],0.005))
  }
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
  pcdata <- pc.res$x
  if (warpgrids == FALSE) {
    if(legend==TRUE){ layout(t(matrix(c(1, 1, 2, 1, 1, 1, 1, 1, 1), 3,3))) }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
         ylab = paste("PC ",axis2))
    if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (length(label!=0)) {
      if(isTRUE(label)){text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
      else{text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.7)) }
    }
    if(!is.null(groups) && legend==TRUE){
      plot.new(); 
      if(is.factor(groups)){legend(0.5,1, legend=unique(groups), pch=19, bty="n", col=unique(groups))
      } else {legend(0.5,1, legend=unique(names(groups)), pch=19, bty="n", col=unique(groups)) }
    }
  }
  shapes <- shape.names <- NULL
  for(i in 1:ncol(pcdata)){
    pcaxis.min <- min(pcdata[, i]) ; pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min ; pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes,pc.min, pc.max)
    shape.names <- c(shape.names,paste("PC",i,"min", sep=""),paste("PC",i,"max", sep=""))
  }
  shapes <- arrayspecs(shapes,p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[,,x])
  names(shapes) <- shape.names
  if (warpgrids == TRUE) {
    if (k == 2) {
      layout(t(matrix(c(2, 1, 4, 1, 1, 1, 1, 1, 3), 3,3)))
    }
    plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21,bg = "black", cex = 2, xlab = paste("PC ",axis1),
         ylab = paste("PC ",axis2))
    if(!is.null(groups)){points(pcdata[, axis1], pcdata[, axis2],pch=21,bg=groups,cex=2)}
    segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
    segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
    if (length(label!=0)) {
      if(isTRUE(label)){text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) }
      else{text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.1)) }
    }
    shape.min <- shapes[[which(names(shapes) == paste("PC",axis1,"min", sep=""))]]
    shape.max <- shapes[[which(names(shapes) == paste("PC",axis1,"max", sep=""))]]
    if (k == 2) {
      arrows(min(pcdata[, axis1]), (0.7 * max(pcdata[,axis2])), min(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      arrows(max(pcdata[, axis1]), (0.7 * min(pcdata[,axis2])), max(pcdata[, axis1]), 0, length = 0.1,lwd = 2)
      tps(ref, shape.min, 20)
      tps(ref, shape.max, 20)
    }
    if(!is.null(groups) && legend==TRUE){
      plot.new(); 
      if(is.factor(groups)){legend(0.5,1, legend=unique(groups), pch=19, bty="n", col=unique(groups))
      } else {legend(0.5,1, legend=unique(names(groups)), pch=19, bty="n", col=unique(groups)) }
    }
    if (k == 3) {
      if (is.null(mesh)==TRUE){
        open3d() ; mfrow3d(1, 2) 
        plot3d(shape.min, type = "s", col = "gray", main = paste("PC ", axis1," negative"),size = 1.25, aspect = FALSE,xlab="",ylab="",zlab="",box=FALSE, axes=FALSE)
        plot3d(shape.max, type = "s", col = "gray", main = paste("PC ", axis1," positive"), size = 1.25, aspect = FALSE,xlab="",ylab="",zlab="",box=FALSE, axes=FALSE)
      }
      if(is.null(mesh)==FALSE){
        open3d() ; mfrow3d(1, 2) 
        cat(paste("\nWarping mesh to negative end of axis ", axis1, "\n", sep=""))
        plotRefToTarget(ref, shape.min, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," negative"))
        next3d()
        cat(paste("\nWarping mesh to positive end of axis ", axis1, "\n", sep=""))
        plotRefToTarget(ref, shape.max, mesh, method = "surface")
        title3d(main=paste("PC ", axis1," positive"))
      }
    }
    layout(1)
  }
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation)
  class(out) = "plotTangentSpace"
  out
}
