# Support functions for gm.prcomp
# AK - 20NOV2017

# Summary
summary.gm.prcomp <- function(x){
  cat(paste("Method: ", attributes(x)$method, "\n", sep=""))
  x$pc.summary
}

# Plot - no warp option, use picknplot for that
### WORK IN PROGRESS - only works for raw pca for now
plot.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, phylo = FALSE, pickNplot = FALSE, ...) {
  pcdata <- x$pc.scores
  plot.args <- list(...)
  plot.args$x <- pcdata[, axis1]  
  plot.args$y <- pcdata[, axis2]
  plot.args$asp <- 1  # Argument to be forced (all others can be handled by the user)
  if(is.null(plot.args$xlab)) plot.args$xlab <- paste("PC ", axis1) 
  if(is.null(plot.args$ylab)) plot.args$ylab <- paste("PC ", axis2)
  if(is.null(plot.args$pch)) plot.args$pch <- 21
  if(is.null(plot.args$bg)) plot.args$bg <- "black"
  
  do.call(plot, args = plot.args)
  segments(0.95*par()$usr[1], 0, 0.95*par()$usr[2], 0, lty = 2, lwd = 1)
  segments(0, 0.95*par()$usr[3], 0, 0.95*par()$usr[4], lty = 2, lwd = 1)

  if(phylo == TRUE){
    if(attributes(x)$method == "Raw data PCA") {
      stop("Raw-data PCA does not allow the projection of the phylogeny.
           Please use phylomorphospace or phylogenetic-PCA in gm.prcomp instead.")
    }
    phy <- attributes(x)$phy
    
    ### Do we need to apply the phylogenetic mean adjustment? And if so, to phylomorphospace only, or also to ppca?
    # pcdata <- pcdata - matrix(x$anc.pcscores[1,], nrow = nrow(pcdata), ncol = ncol(pcdata), byrow = T)
    
    
    
  }  
}
