# Function for ace of GM data
# follows fastAnc in phytools

shape.ace <- function(x, phy){
  N <- length(phy$tip.label)
  Nnode <- phy$Nnode
  anc.states<-NULL   
  for (i in 1:ncol(x)){
    x1 <- x[,i]
    tmp <- vector()
    for (j in 1:Nnode + N) {
      a <- multi2di(root(phy, node = j))
      tmp[j - N] <- ace(x1, a, method = "pic")$ace[1]
      }
    anc.states<-cbind(anc.states,tmp)   
    }
  colnames(anc.states) <- NULL
  row.names(anc.states) <- 1:length(tmp)
  return(anc.states)
}
