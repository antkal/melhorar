library(geomorph)

rm(list=ls())
source("gm.prcomp.R")
source("support.gm.prcomp.R")

# Example data (with phylo)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)


# Testing gm.prcomp
test.phy <- gm.prcomp(A=Y.gpa$coords, phy=plethspecies$phy, phylo.pca = TRUE)
plot(test.phy$pc.scores, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))

library(phytools)
test.phy2 <- phyl.pca(tree=plethspecies$phy, Y=two.d.array(Y.gpa$coords))
plot(test.phy2$S, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))
# not same. What are we doing differently?

layout(matrix(1:2, nrow=2))
plot(test.phy$pc.scores, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))
plot(test.phy2$S, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))

d1 <- dist(test.phy$pc.scores)
d2 <- dist(test.phy2$S)
plot(d1, d2)
