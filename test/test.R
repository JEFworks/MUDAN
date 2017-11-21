library(MUDAN)

## test functions using simulated data
simulate.data <- function(G=5, N=100, M=1000, initmean=0, initvar=10, upreg=10, upregvar=10, ng=100, seed=0, plot=TRUE) {
  set.seed(seed)
  mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
  rownames(mat) <- paste0('gene', 1:M)
  colnames(mat) <- paste0('cell', 1:(N*G))
  group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
  names(group) <- colnames(mat)

  diff <- lapply(1:G, function(x) {
    diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
    mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
    return(diff)
  })
  names(diff) <- paste0('group', 1:G)

  diff2 <- lapply(2:(G-1), function(x) {
    y <- x+G
    diff <- rownames(mat)[(((y-1)*ng)+1):(((y-1)*ng)+ng)]
    mat[diff, group %in% paste0("group", 1:x)] <<- mat[diff, group %in% paste0("group", 1:x)] + rnorm(ng, upreg, upregvar)
    return(diff)
  })

  mat[mat<0] <- 0
  mat <- round(mat)

  if(plot) {
    heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)
  }

  return(list(mat=mat, group=group))
}

## simulate data
data <- simulate.data()
group <- data$group
set.seed(0)
## split into two datasets
d1 <- c(sample(1:100, 90), sample(101:200, 90), sample(201:300, 50), sample(301:400, 50), sample(401:500, 5))
d2 <- setdiff(1:length(group), d1)
d1cells <- names(group)[d1]
d2cells <- names(group)[d2]
table(group[d1cells])
table(group[d2cells])

## first group has rare purple cell type that is difficult to distinguish
myMudanObject1 <- Mudan$new("d1", data$mat[, d1cells])
myMudanObject1$libSizeNormalize()
myMudanObject1$varianceNormalize(plot=FALSE)
myMudanObject1$dimensionalityReduction(maxit=1000)
myMudanObject1$getStandardEmbedding(groups=group[d1cells], plot=TRUE)
myMudanObject1$communityDetection(reductionType='pcs', communityName="Infomap", communityMethod=igraph::cluster_infomap, k=10)
myMudanObject1$modelCommunity(communityName="Infomap")
myMudanObject1$getMudanEmbedding()
myMudanObject1$plot(reductionType='pcs', communityName="Infomap", embeddingType="MUDAN")
myMudanObject1$plot(groups=group[d1cells], embeddingType="MUDAN")

## second group has much easier time distinguishing purple cell type
myMudanObject2 <- Mudan$new("d2", data$mat[, d2cells])
myMudanObject2$libSizeNormalize()
myMudanObject2$varianceNormalize(plot=FALSE)
myMudanObject2$dimensionalityReduction(maxit=1000)
myMudanObject2$getStandardEmbedding(groups=group[d2cells], plot=TRUE)
myMudanObject2$communityDetection(reductionType='pcs', communityName="Infomap", communityMethod=igraph::cluster_infomap, k=10)
myMudanObject2$modelCommunity(communityName="Infomap")
myMudanObject2$getMudanEmbedding()
myMudanObject2$plot(reductionType='pcs', communityName="Infomap", embeddingType="MUDAN")
myMudanObject2$plot(groups=group[d2cells], embeddingType="MUDAN")

## learn from myMudanObject2
reduction <- predict(myMudanObject2$model, data.frame(t(as.matrix(myMudanObject1$mat))))$x
emb <- Rtsne::Rtsne(lds12, is_distance=FALSE, perplexity=30, verbose=TRUE, num_threads=2)$Y
rownames(emb) <- rownames(reduction)
plotEmbedding(emb, groups=group[rownames(emb)]) ## now purple comes out

## but myMudanObject1 has its own advantages, so why not leverage both
## just combine embeddings
lds11 <- predict(myMudanObject1$model, data.frame(t(as.matrix(myMudanObject1$mat))))$x
lds12 <- predict(myMudanObject2$model, data.frame(t(as.matrix(myMudanObject1$mat))))$x
lds21 <- predict(myMudanObject1$model, data.frame(t(as.matrix(myMudanObject2$mat))))$x
lds22 <- predict(myMudanObject2$model, data.frame(t(as.matrix(myMudanObject2$mat))))$x
reduction <- rbind(cbind(lds11, lds12), cbind(lds21, lds22))
emb <- Rtsne::Rtsne(reduction, is_distance=FALSE, perplexity=30, verbose=TRUE, num_threads=2)$Y
rownames(emb) <- rownames(reduction)
plotEmbedding(emb, groups=group[rownames(emb)])
