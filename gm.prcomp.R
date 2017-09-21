library(geomorph)

# Example data (with phylo)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)

source("geomorph.support.code.r")

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
  p <- dim(A)[1]; k <- dim(A)[2]; n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)

  if(!is.null(phy)){
    if (!inherits(phy, "phylo"))
      stop("Tree must be of class 'phylo.'")
    if (!is.binary.tree(phy)) 
      stop("Tree is not fully bifurcating (consider 'multi2di' in ape.")
    N <- length(phy$tip.label); Nnode <- phy$Nnode
    if(N!=n)
      stop("Number of taxa in data matrix and tree are not equal.")
    if(is.null(rownames(x))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
    } else x <- x[phy$tip.label, ]  
    
    # PhyloPCA (follows phylo.integration code)
    phy.parts <- phylo.mat(x,phy)
    invC <- phy.parts$invC; D.mat <- phy.parts$D.mat
    
    one <- matrix(1,nrow(x)); I <- diag(1,nrow(x)) 
    Ptrans <- D.mat%*%(I-one%*%crossprod(one,invC)/sum(invC))
    x <- Ptrans%*%x
    pca <- La.svd(x)
    PCscores <- x %*% t(pca$vt)
  }
  
  if(is.null(tol)){
    d <- prcomp(x)$sdev^2
    cd <-cumsum(d)/sum(d)
    cd <- length(which(cd < 1)) 
    if(length(cd) < length(d)) cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }

  
  
  
  pc.res <- prcomp(x, center = center, scale. = scale., retx = retx, tol = tol)
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
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata, pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation)
  class(out) = "gm.prcomp"
  out
}
