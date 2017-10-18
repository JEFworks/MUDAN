############ Developing faster KNN cluster membership using subsampling
library(MUDAN)
source('../R/mudan.R')

library(Matrix)
data(pbmcA)
cd <- cleanCounts(pbmcA)
cd.norm <- normalizeCounts(cd)
cd.varnorm <- log10(normalizeVariance(cd.norm)+1)
pcs <- getPcs(cd.varnorm, maxit=1000)
dim(pcs)

############ Function development

## ## random subsample
## subsample <- sample(colnames(pcs), 1000)
## pcs.sub <- pcs[, subsample]
## com.sub <- getKnnMembership(pcs.sub, k=10, method=igraph::cluster_infomap)

## ## check what it looks like
## pcs.emb <- Rtsne::Rtsne(t(pcs), is_distance=FALSE, perplexity=30, num_threads=10)$Y
## rownames(pcs.emb) <- colnames(pcs)
## plotEmbedding(pcs.emb,
##               groups=com.sub,
##               mark.clusters = FALSE,
##               show.legend=TRUE,
##               legend.x = 'bottomleft',
##               alpha=0.1)

## ## for unannotated cells, look for nearest neighbor in annotated cells
## data <- pcs[, subsample]
## query <- pcs[, setdiff(colnames(pcs), subsample)]
## knn <- RANN::nn2(t(data), t(query), k=30)[[1]]
## rownames(knn) <- colnames(query)

## com.nonsub <- unlist(apply(knn, 1, function(x) {
##     ## nearest neighbors in data
##     nn <- colnames(data)[x]
##     ## look at their cell type annotations
##     nn.com <- com.sub[nn]
##     ## get most frequent annotation
##     return(names(sort(table(nn.com), decreasing=TRUE)[1]))
## }))

################# Cleaned up
com.approx <- getApproxKnnMembership(pcs, 10, 10, nsubsample=1000, method=igraph::cluster_infomap)
par(mfrow=c(1,2))
plotEmbedding(pcs.emb,
              groups=com.approx,
              mark.clusters = FALSE,
              show.legend=TRUE,
              legend.x = 'bottomleft',
              alpha=0.2)

## compared to running on full dataset
com.all <- getKnnMembership(pcs, k=15, method=igraph::cluster_infomap)
plotEmbedding(pcs.emb,
              groups=com.all,
              mark.clusters = FALSE,
              show.legend=TRUE,
              legend.x = 'bottomleft',
              alpha=0.2)

table(com.all, c(com.sub, com.nonsub)[names(com.all)])
