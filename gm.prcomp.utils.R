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
}
