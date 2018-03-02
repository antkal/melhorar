# Support functions for gm.prcomp
# AK - 20NOV2017
source("plotRefToTarget.r")

# Summary
summary.gm.prcomp <- function(x){
  cat(paste("Method: ", attributes(x)$method, "\n", sep=""))
  x$pc.summary
}

# Plot - no warp option, use picknplot for that
plot.gm.prcomp <- function(X, axis1 = 1, axis2 = 2, phylo = FALSE, 
                           phylo.par = list(edge.color = "black", edge.width = 1, edge.lty = 1,
                           node.bg = "black", node.pch = 21, node.cex = 1), ...) {
  pcdata <- X$pc.scores
  plot.args <- list(...)
  plot.args$x <- pcdata[, axis1]  
  plot.args$y <- pcdata[, axis2]
  plot.args$asp <- 1  
  if(is.null(plot.args$xlab)) plot.args$xlab <- paste("PC ", axis1) 
  if(is.null(plot.args$ylab)) plot.args$ylab <- paste("PC ", axis2)

  if(phylo == FALSE){
    do.call(plot, args = plot.args)
  } else {
    if(attributes(X)$method == "Raw data PCA") {
      stop("Raw-data PCA does not allow the projection of a phylogeny.
           Please use phylomorphospace or phylogenetic-PCA.")
    }
    pcdata <- rbind(pcdata, X$anc.pcscores)
    phy <- attributes(X)$phy
    ### Do we need to apply the phylogenetic mean adjustment? 
    # And if so, to phylomorphospace only, or also to ppca?
    # pcdata <- pcdata - matrix(x$anc.pcscores[1,], nrow = nrow(pcdata), ncol = ncol(pcdata), byrow = T)
    plot.args$type <- "n"
    do.call(plot, args = plot.args)
    for (i in 1:nrow(phy$edge)) {
      lines(pcdata[(phy$edge[i,]), axis1], pcdata[(phy$edge[i,]), axis2], type="l",
            col = phylo.par$edge.color, lwd = phylo.par$edge.width, lty = phylo.par$edge.lty)
    }
    plot.args$type <- "p"
    do.call(points, args = plot.args)
    points(pcdata[(length(phy$tip)+1):nrow(pcdata), axis1], pcdata[(length(phy$tip)+1):nrow(pcdata), axis2],
           pch = phylo.par$node.pch, cex = phylo.par$node.cex, bg = phylo.par$node.bg)
  }
  continue <- readline("Do you want to visualize shape variation in morphospace (y/n)? ")
  while(continue == "y"){
    cat("Pick a point in shape space", "\n")
    picked.pts <- unlist(locator(n = 1, type = "p", pch = 20, col = "red", cex = 1))
    cat("Picked point coordinates are:", "\n")
    cat(picked.pts, "\n") 
    preds <- shape.predictor(attributes(X)$Adata, x = X$pc.scores[,c(axis1, axis2)], Intercept = FALSE, pred1 = picked.pts)
    if (dim(attributes(X)$Adata)[2]==2) {
      plot.args$M1 <- cbind(mshape(attributes(X)$Adata), 0)
      plot.args$M2 <- cbind(preds$pred1, 0)
      class(plot.args$M2) <- "predshape.k2"
      view3d(phi = 0, fov = 30, interactive = FALSE) 
      do.call(plotRefToTarget, plot.args)
    }
    if (dim(attributes(X)$Adata)[2]==3){
      plot.args$M1 <- mshape(attributes(X)$Adata)
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
