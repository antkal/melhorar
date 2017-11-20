# gm.prcomp - examples and tests
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
gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4)))

plot(pleth.raw)
plot(pleth.raw, pch=22, cex = 1.5, xlab = "Maria", asp = NULL, bg = gps) # Modify options
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




