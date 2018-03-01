#' Principal components analysis of shape data
#'
#' Function performs raw or weighted PCA on superimposed shape coordinates  
#'
#' The function performs a principal components analysis of shape variation, with the possibility
#' of weighing  a phylogenetic tree, or a variance-covariance matrix, to incorporate the expected
#'   
#' 
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of aligned specimens 
#' @param warpgrids A logical value indicating whether deformation grids for shapes along X-axis should be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation along X-axis (when {warpgrids=TRUE})
#' as described in \code{\link{plotRefToTarget}}.
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param label An optional vector indicating labels for each specimen are to be displayed 
#' (or if TRUE, numerical addresses are given)
#' @param groups An optional factor vector specifying group identity for each specimen (see example)
#' @param legend A logical value for whether to add a legend to the plot (only when groups are assigned).
#' @param ... Arguments passed on to \code{\link{prcomp}}.  By default, \code{\link{plotTangentSpace}}
#' will attempt to remove redundant axes (eigen values effectively 0).  To override this, adjust the 
#' argument, tol, from \code{\link{prcomp}}.
#' @return If user assigns function to object, returned is a list of the following components:
#' \item{pc.summary}{A table summarizing the percent variation explained by each pc axis, equivalent to summary of \code{\link{prcomp}}.}
#' \item{pc.scores}{The set of principal component scores for all specimens.}
#' \item{pc.shapes}{A list with the shape coordinates of the extreme ends of all PC axes, e.g. $PC1min}
#' \item{sdev}{The standard deviations of the principal components (i.e., the square roots of the eigenvalues of the 
#' covariance/correlation matrix, as per \code{\link{prcomp}}.}
#' \item{rotation}{The matrix of variable loadings, as per \code{\link{prcomp}}.}
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' 
#' gp <- interaction(plethodon$species, plethodon$site) # group must be a factor
#' plotTangentSpace(Y.gpa$coords, groups = gp) 
#' 
#' ## To save and use output
#' PCA <- plotTangentSpace(Y.gpa$coords, groups = gp, legend=TRUE) 
#' summary(PCA)
#' PCA$pc.shapes
#' PCA$rotation
#' 
#' ##To change colors of groups
#' col.gp <- rainbow(length(levels(gp))) 
#'    names(col.gp) <- levels(gp)
#' col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must NOT be a factor
#' plotTangentSpace(Y.gpa$coords, groups = col.gp)
#' 
#' ## To plot residual shapes from an allometry regression (note: must add mean back in!) 
#' plotTangentSpace(arrayspecs(resid(lm(two.d.array(Y.gpa$coords)~Y.gpa$Csize))+
#'          predict(lm(two.d.array(Y.gpa$coords)~1)),12,2))

source("shape.ace.R")
source("cov.mat.R")
source("geomorph.support.code.r")

# Only A: normal, raw PCA, accepts pca arguments through ...
# A + phy: GMphylomorphospace
# A + phy + phylo.pca = T: phyloPCA
# A + COV: other weighed PCA (with catch for phy&cov)

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
  
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata[1:n, ], pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation)
  
  if(meth == "Phylomorphospace") {
    out$anc.states <- anc
    out$anc.pcscores <- pcdata[(N+1):nrow(pcdata),]
  }
  
  if(meth == "Phylogenetic PCA"){
    out$anc.states <- anc
    out$anc.pcscores <- anc%*%pc.res$rotation
  }
  
  class(out) = "gm.prcomp"
  attributes(out)$method <- meth
  attributes(out)$phy <- phy
  attributes(out)$Adata <- A
  return(out)
}
