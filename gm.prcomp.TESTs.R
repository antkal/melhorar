# gm.prcomp - examples and tests
rm(list=ls())
source("gm.prcomp.R")
source("support.gm.prcomp.R")

# Example data (with phylo)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)

# PCAs
pleth.raw <- gm.prcomp(Y.gpa$coords)
pleth.phylomorpho <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
pleth.ppca <- gm.prcomp(Y.gpa$coords, plethspecies$phy, phylo.pca = T)

# Summaries
summary(pleth.raw)
summary(pleth.phylomorpho)
summary(pleth.ppca) 

# Plot of raw PCA data
gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups

plot(pleth.raw) #Note that it prompts a question for visualizing shapes that calls incorporated picknplot.shape code
par(mgp = c(2.5, 0.5, 0))
plot(pleth.raw, pch=22, cex = 1.5, xlab = "My PCA - axis 1", asp = NULL, bg = gps,
     font.lab = 2, cex.lab = 2) # Modify options - asp is forced to asp = 1
# Add things as desired using standard R plotting
segments(0.95*par()$usr[1], 0, 0.95*par()$usr[2], 0, lty = 2, lwd = 1)
segments(0, 0.95*par()$usr[3], 0, 0.95*par()$usr[4], lty = 2, lwd = 1)
legend("topright", pch=22, pt.bg = unique(gps), legend = levels(gps), cex = 2) # Modify options

# Plot of phylomorphospace and phyloPCA scores, without tree
plot(pleth.phylomorpho)
text(pleth.phylomorpho$pc.scores, labels = labels(pleth.phylomorpho$pc.scores)[[1]],
     adj = c(-0.1, -0.1)) # The usual issues with labeling, fixed by hand e.g. like this:

plot(pleth.phylomorpho, xlim=c(-0.035, 0.04))
text(pleth.phylomorpho$pc.scores, labels = labels(pleth.phylomorpho$pc.scores)[[1]],
     adj = c(-0.1, -0.1)) 

plot(pleth.ppca, xlim=c(-0.01, 0.015), cex=1.5, pch=23, bg=ifelse(gps=="gp1", "red", "blue"))
text(pleth.ppca$pc.scores, labels = labels(pleth.ppca$pc.scores)[[1]],
     adj = c(-0.1, -0.1)) 
legend("bottomleft", pch=23, pt.bg = c("red", "blue"), legend = levels(gps), cex = 1.5, bty = "n")

# Plots WITH trees
plot(pleth.raw, phylo = T)  # Gives error (as it should)
plot(pleth.phylomorpho, phylo = T)

# Modify parameters - arguments passed directly modify tip plotting
plot(pleth.phylomorpho, phylo = T, cex = 2, pch = 22, bg = gps, lwd = 2)
plot(pleth.phylomorpho, phylo = T, cex = 2, pch = 22, bg = gps, 
     phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
                      node.cex = 0)) # Supresses plotting of nodes

# Plot and add things by hand (base R plotting)
par(mgp = c(2, 0.5, 0))
plot(pleth.ppca, phylo = T, cex = 1.5, pch = 22, bg = gps, cex.lab = 2, font.lab = 2, xlim = c(-0.007, 0.017),
     phylo.par = list(edge.color = "grey", edge.width = 2,
                      node.bg = "black", node.pch = 22, node.cex = 0.5)
     )
text(pleth.ppca$pc.scores, labels = labels(pleth.ppca$pc.scores)[[1]],
     adj = c(-0.1, -0.1), font = 4) 
text(pleth.ppca$anc.pcscores, labels = labels(pleth.ppca$anc.pcscores)[[1]],
  adj = c(-0.1, -0.1), font = 2) 

# Does it make sense to allow plotting ppcas with the tree projected??? 
# Phylogeny is "standardized for" in this version
