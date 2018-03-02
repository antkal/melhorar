# gm.prcomp

#' Print/Summary function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
print.gm.prcomp <- function (x, ...) {
  cat(paste("Method: ", attributes(x)$method, "\n", sep=""))
  print(x$pc.summary)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
summary.gm.prcomp <- function (object, ...) {
  print.gm.prcomp(object, ...)
}

#' Plot Function for geomorph
#' 
#' To visualize shape variation across PC axes, use \code\{\link{plotRefTotarget}} with the
#' $pc.shapes component of your gm.prcomp object.
#' 
#' @param x An object of class \code{\link{gm.prcomp}}
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param phylo A logical value indicating whether the phylogeny should be projected to PC space
#' @param phylo,par A list of plotting parameters for the phylogeny edges (edge.color, edge.width, edge.lty)
#' and nodes (node.bg, node.pch, node.cex)
#' @param ... other arguments passed to plot
#' @return An object of class "plot.gm.prcomp" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
#' @keywords visualization
#' @seealso  \code{\link{plotRefToTarget}}


plot.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, phylo = FALSE, 
                           phylo.par = list(edge.color = "black", edge.width = 1, edge.lty = 1,
                           node.bg = "black", node.pch = 21, node.cex = 1), ...) {
  pcdata <- x$pc.scores[, c(axis1, axis2)]

  if(phylo == FALSE){
    plot.new()
    plot.window(1.05*range(pcdata[,1]), 1.05*range(pcdata[,2]), asp=1,...)
    plot.xy(xy.coords(pcdata), type="p",...)
  } else {
    if(attributes(x)$method == "Raw data PCA") {
      stop("Raw-data PCA does not allow the projection of a phylogeny.
           Please use phylomorphospace or phylogenetic-PCA.")
    }
    phy <- attributes(x)$phy
    pcdata <- rbind(pcdata, x$anc.pcscores[,c(axis1, axis2)])
    plot.new()
    plot.window(1.05*range(pcdata[,1]), 1.05*range(pcdata[,2]), asp=1,...)
    for (i in 1:nrow(phy$edge)) {
      dt.xy <- xy.coords(pcdata[(phy$edge[i,]),])
      plot.xy(dt.xy, type="l", col = phylo.par$edge.color, 
              lwd = phylo.par$edge.width, lty = phylo.par$edge.lty,...)
    }
    plot.xy(xy.coords(pcdata[1:length(phy$tip),]), type="p",...)
    plot.xy(xy.coords(pcdata[(length(phy$tip)+1):nrow(pcdata),]), type="p",
            pch = phylo.par$node.pch, cex = phylo.par$node.cex, bg = phylo.par$node.bg,...)
  }
  out <- list(points = x$pc.data[,1:2], pc.data = x$pc.data)
  class(out) <- "plot.gm.prcomp"
  invisible(out)
  
}
