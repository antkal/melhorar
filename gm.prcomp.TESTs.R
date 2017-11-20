# gm.prcomp - examples and tests
source("gm.prcomp.R")

# Example data (with phylo)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)

pleth.raw <- gm.prcomp(Y.gpa$coords)
pleth.phylomorpho <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
pleth.ppca <- gm.prcomp(Y.gpa$coords, plethspecies$phy, phylo.pca = T)

summary(pleth.raw)
summary(pleth.phylomorpho)
summary(pleth.ppca)

# Scores
x.raw <- pleth.raw$pc.scores[plethspecies$phy$tip.label,]
x.pms <- pleth.phylomorpho$pc.scores
x.ppca <- pleth.ppca$pc.scores

d.raw <- dist(x.raw)
d.pms <- dist(x.pms)
d.ppca <- dist(x.ppca)

plot(d.raw, d.pms) # Straight line, cor = 1, OK :)



