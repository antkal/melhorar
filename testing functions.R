# Testing gm.prcomp
test <- gm.prcomp(Y.gpa$coords)
plot(test$pc.scores, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))

test.phy <- gm.prcomp(A=Y.gpa$coords, phy=plethspecies$phy, phylo.pca = TRUE)
plot(test.phy$pc.scores, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))

library(phytools)
test.phy2 <- phyl.pca(tree=plethspecies$phy, Y=two.d.array(Y.gpa$coords))
plot(test.phy2$S, pch=19, col=rainbow(dim(Y.gpa$coords)[3]))
# not same. What are we doing differently?


