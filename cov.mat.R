# Function for retrieving pieces for PCA weighting from a VCV (not from phylo)
# Follows phylo.mat
# AK - 20NOV2017

cov.mat <- function(x, COV) {
  if(is.null(dimnames(COV))){
    warning("COV matrix does not include dimnames. Assuming same order of observations as x")
  } else {
    C <- C[rownames(x),rownames(x)]
  }
  invC <- fast.solve(C)
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect <- eigC$vectors[,1:(length(lambda))]
  D.mat <- fast.solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect))
  rownames(D.mat) <- colnames(D.mat) <- colnames(C)
  rownames(invC) <- colnames(invC) <- colnames(C)
  list(invC = invC, D.mat = D.mat, C = C)
}
