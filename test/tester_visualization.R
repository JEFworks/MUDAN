#################### Compare MUDAN visualization to conventional PCA-based tSNE
library(MUDAN)
source('../R/mudan.R')

library(Matrix)
data(pbmcA)
cd <- cleanCounts(pbmcA)
cd.norm <- normalizeCounts(cd)
cd.varnorm <- log10(normalizeVariance(cd.norm)+1)

pcs <- getPcs(cd.varnorm, maxit=1000)
## try list of ks
ks <- c(10,30,50,100)
coms <- lapply(ks, function(k) {
    com <- getKnnMembership(pcs, k=k, method=igraph::cluster_infomap)
    min.group.size <- 30
    bad.groups <- names(which(table(com) < min.group.size))
    com[com %in% bad.groups] <- NA
    return(com)
})

## mudan embedding
vargenes <- getVariableGenes(cd.varnorm, 1000)
models <- lapply(coms, function(com) {
    model <- modelLda(mat=cd.varnorm[vargenes,], com=com)
})
mudan.emb <- tsneLda(mat=cd.varnorm, model=models, com=com, perplexity=30, plot=FALSE, details=TRUE)
## run new knn
mudan.com <- getKnnMembership(t(mudan.emb$reduction), k=30, method=igraph::cluster_walktrap)

## get PCA based embedding
pcs.emb <- Rtsne::Rtsne(t(pcs), is_distance=FALSE, perplexity=30, num_threads=10)$Y
rownames(pcs.emb) <- colnames(pcs)

## plot to compare
png(filename = "../docs/tester_visualization.png", width = 600, height = 300)
par(mfrow=c(1,2))
plotEmbedding(pcs.emb,
              groups=mudan.com,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='PCA-based tSNE')
plotEmbedding(mudan.emb$emb,
              groups=mudan.com,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='MUDAN')
dev.off()

#################### Try different community detection
mudan.com1 <- getKnnMembership(t(mudan.emb$reduction), k=30, method=igraph::cluster_walktrap)
mudan.com2 <- getKnnMembership(t(mudan.emb$reduction), k=30, method=igraph::cluster_infomap)
mudan.com3 <- getKnnMembership(t(mudan.emb$reduction), k=30, method=igraph::cluster_louvain)
mudan.com4 <- getKnnMembership(t(mudan.emb$reduction), k=30, method=igraph::cluster_fast_greedy)
png(filename = "../docs/tester_coms.png", width = 600, height = 600)
par(mfrow=c(2,2))
plotEmbedding(mudan.emb$emb,
              groups=mudan.com1,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='Walktrap Community Detection')
plotEmbedding(mudan.emb$emb,
              groups=mudan.com2,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='Infomap Community Detection')
plotEmbedding(mudan.emb$emb,
              groups=mudan.com3,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='Louvain Community Detection')
plotEmbedding(mudan.emb$emb,
              groups=mudan.com4,
              mark.clusters = FALSE,
              show.legend=FALSE,
              legend.x = 'topleft',
              alpha=0.1,
              main='Fast Greedy Community Detection')
dev.off()

## diff genes
diffGenes <- getDifferentialGenes(cd.norm, mudan.com)
## markers
markers <- unlist(lapply(diffGenes$info, function(x) {
    rownames(x)[x$marker.auc>0.8]
}))
markers

## plot resulting embedding
heatmap(cd.norm[markers,names(sort(mudan.com))], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", labCol=FALSE, labRow=FALSE, main='\nIdentified markers: AUC > 0.9')

## gene set enrichment of diff genes in group 1
