library(geomorph)

source("geomorph.support.code.r")
source("shape.ace.R")
source("cov.mat.R")

# Only A: normal, raw PCA, accepts pca arguments through ...
# A + phy: GMphylomorphospace
# A + phy + phylo.pca = T: phyloPCA
# A + COV: other weighed PCA (with catch for phy&cov not implemented)

gm.prcomp <- function (A, phy = NULL, phylo.pca = FALSE, COV = NULL, ...){
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
  p <- dim(A)[1]; k <- dim(A)[2]; n <- dim(A)[3]
  ref <- mshape(A)
  x <- scale(two.d.array(A), center = center, scale = scale.)

  if(is.null(phy) & phylo.pca == T){
    stop("To perform phylogenetic pca, please provide a phylogeny.")
  }
  
  if(!is.null(phy)){
    if (!inherits(phy, "phylo"))
      stop("Tree must be of class 'phylo.'")
    if (!is.binary.tree(phy)) 
      stop("Tree is not fully bifurcating (consider 'multi2di' in ape).")
    if(!is.null(COV)){
      stop("Method not implemented for weighting with BOTH a phylogeny and another covariance matrix.")
    }
    N <- length(phy$tip.label); Nnode <- phy$Nnode
    if(N!=n)
      stop("Number of taxa in data matrix and tree are not equal.")
    if(is.null(rownames(x))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
    } else x <- x[phy$tip.label, ]
    anc <- shape.ace(x, phy)
    
    if(phylo.pca == T){
      # Phylogenetic transformation (follows phylo.integration code)
      phy.parts <- phylo.mat(x, phy)
      invC <- phy.parts$invC; D.mat <- phy.parts$D.mat
      one <- matrix(1, nrow(x)); I <- diag(1, nrow(x)) 
      Ptrans <- D.mat%*%(I-one%*%crossprod(one, invC)/sum(invC))
      x <- Ptrans%*%x
    } else {
      x <- rbind(x, anc)
    }
  } else {
    if(!is.null(COV)){
      # Transformation of data using COV 
      # (as per phylo.integration, but using a COV matrix instead of phylo)
      cov.parts <- cov.mat(x, COV)
      invC <- cov.parts$invC; D.mat <- cov.parts$D.mat
      one <- matrix(1, nrow(x)); I <- diag(1, nrow(x)) 
      Ptrans <- D.mat%*%(I-one%*%crossprod(one, invC)/sum(invC))
      x <- Ptrans%*%x
    }
  }
  
  if(is.null(tol)){
    d <- prcomp(x)$sdev^2
    cd <-cumsum(d)/sum(d)
    cd <- length(which(cd < 1)) 
    if(length(cd) < length(d)) cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }

  pc.res <- prcomp(x, retx = retx, tol = tol)
  pcdata <- pc.res$x
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
  
  # OUTPUT
  if(is.null(phy) & is.null(COV)) meth <- "Raw data PCA"
  if(!is.null(phy) & phylo.pca == FALSE) meth <- "Phylomorphospace"
  if(!is.null(phy) & phylo.pca == TRUE) meth <- "Phylogenetic PCA"
  if(!is.null(COV)) meth <- "COV-weighted PCA"
  
  if(!is.null(phy) | !is.null(COV)) pcscores <- pcdata[1:n, ]

  out <- list(pc.summary = summary(pc.res), pc.scores = pcscores, pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation)
  if(!is.null(phy)) {
    out$anc.states <- anc
    out$anc.pcscores <- pcdata[(N+1):nrow(pcdata),]
  }
  class(out) = "gm.prcomp"
  attributes(out)$method <- meth
  attributes(out)$phy <- phy
  return(out)
}
