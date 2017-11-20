# Support functions for gm.prcomp
# AK - 20NOV2017

# Summary
summary.gm.prcomp <- function(x){
  cat(paste("Method: ", x$method, "\n", sep=""))
  x$pc.summary
}

# Plot - no warp option, use picknplot for that
### WORK IN PROGRESS
plot.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, phylo = TRUE, 
                           groups = NULL, legend = FALSE, ...){
  m <- x$method
  #if(m == "Raw data PCA")
  
  pcdata <- x$pc.scores
  if(legend == TRUE) { 
    layout(t(matrix(c(1, 1, 2, 1, 1, 1, 1, 1, 1), 3,3))) 
    }
  plot(pcdata[, axis1], pcdata[, axis2], asp = 1, pch = 21, bg = "black", cex = 2, xlab = paste("PC ", axis1),
       ylab = paste("PC ", axis2))
  if(!is.null(groups)) {
    points(pcdata[, axis1], pcdata[, axis2], pch=21, bg=groups, cex=2)
    }
  segments(min(pcdata[, axis1]), 0, max(pcdata[, axis1]), 0, lty = 2, lwd = 1)
  segments(0, min(pcdata[, axis2]), 0, max(pcdata[, axis2]), lty = 2, lwd = 1)
  if (length(label!=0)) {
    if(isTRUE(label)) { 
      text(pcdata[, axis1], pcdata[, axis2], seq(1, n), adj = c(-0.7, -0.7)) 
    } else { 
        text(pcdata[, axis1], pcdata[, axis2], label, adj = c(-0.1, -0.7)) 
      }
  }
  if(!is.null(groups) && legend==TRUE){
    plot.new(); 
    if(is.factor(groups)){legend(0.5,1, legend=unique(groups), pch=19, bty="n", col=unique(groups))
    } else {legend(0.5,1, legend=unique(names(groups)), pch=19, bty="n", col=unique(groups)) }
  }
  
}
