library(MUDAN)
data("pbmcA")

myMudanObject1 <- Mudan$new("pbmcA", pbmcA, ncores=4)
myMudanObject1$libSizeNormalize()
myMudanObject1$varianceNormalize(plot=FALSE)
myMudanObject1$dimensionalityReduction(nGenes = 1000, nPcs = 30, maxit=1000)
myMudanObject1$getStandardEmbedding(plot=FALSE)

## community detection: overcluster
myMudanObject1$communityDetection(reductionType='pcs', communityName="Infomap", communityMethod=igraph::cluster_infomap, k=10)
myMudanObject1$plot(reductionType='pcs', communityName="Infomap", embeddingType="PCA", main='pbmcA PCA', show.legend=TRUE, mark.clusters=TRUE)
cols <- na.omit(myMudanObject1$com[['pcs']][['Infomap']])
cm <- t(myMudanObject1$cd[, names(cols)])
cols2 <- getStableClusters(cm, cols)

## if super overclustered may need to recur?
cols3 <- getStableClusters(cm, cols2)
cols4 <- getStableClusters(cm, cols3)
cols5 <- getStableClusters(cm, cols4)
