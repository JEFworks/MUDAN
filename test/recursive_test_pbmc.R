library(Matrix)
data(pbmcA)
cd <- normalizeCounts(pbmcA)
class(cd)
cd0 <- cd
matnorm <- normalizeVariance(cd, plot=TRUE, details=TRUE, alpha=0.05)
names(matnorm)
mat <- log10(matnorm$mat+1)[matnorm$ods,]
dim(mat)
mat0 <- mat

nPcs = 30
pcs <- fastPca(t(mat), nPcs)
m <- pcs$l
colnames(m) <- paste0('PC', seq_len(nPcs))
##variance explained
#vars_transformed <- apply(m, 2, var)
#ve <- sort(vars_transformed, decreasing=TRUE)
#plot(ve, type="l")
#abline(h=mean(ve), col='red')
#pcs.fin <- m[, names(which(ve>mean(ve)))]
pcs.fin <- m
#pcs.fin <- m[, 1:max(which(cumsum(ve)<5))] ## require 5% of variance be explained
#dim(pcs.fin)

emb1 <- Rtsne::Rtsne(as.matrix(pcs.fin), is_distance=FALSE, perplexity=30, verbose=TRUE, num_threads=4)$Y
rownames(emb1) <- colnames(mat)

#com0 <- getComMembership(t(as.matrix(pcs.fin)), k=30, method=igraph::cluster_walktrap)
com0 <- getComMembership(t(as.matrix(pcs.fin)), k=100, method=igraph::cluster_walktrap)
plotEmbedding(emb1, groups=com0, show.legend=TRUE, main='Overclustered')
com0.stable <- getStableClusters(cm=t(cd0), cols=com0)
plotEmbedding(emb1, groups=com0.stable, show.legend=TRUE, main='Stable')
#com0.stable2 <- getStableClusters(cm=t(cd0), cols=com0.stable)
#plotEmbedding(emb1, groups=com0.stable2, show.legend=TRUE, main='Stable')

ds <- getDifferentialGenes(cm=t(cd0), cols=com0.stable)
diff.genes <- lapply(ds, function(x) rownames(x)[x$Z>1.96])
lapply(diff.genes, na.omit)
na.omit(diff.genes[[7]])
na.omit(diff.genes[[10]])
na.omit(diff.genes[[11]])

model <- modelLda(mat0, com0.stable)
reduction <- predict(model, data.frame(t(as.matrix(mat))))$x
emb2 <- Rtsne::Rtsne(reduction, is_distance=FALSE, perplexity=30, verbose=TRUE, num_threads=4)$Y
rownames(emb2) <- colnames(mat)
par(mfrow=c(1,2))
plotEmbedding(emb1, groups=com0.stable, show.legend=TRUE, main='PCA')
plotEmbedding(emb2, groups=com0.stable, show.legend=TRUE, main='LDA')

## plot some genes
gs <- c('CD3D', 'CD4', 'CD8A', 'NKG7', 'FCER1A', 'FCGR3A', 'S100A8',
        'CCR10', 'TNFRSF18', 'FOXP3', 'ID3', 'GZMB', 'MS4A1', 'CD19',
        'GP1BA', 'GZMK', 'GZMA')
gs <- intersect(gs, rownames(cd0))
length(gs)
par(mfrow=c(4,4))
invisible(lapply(gs, function(g) {
  plotEmbedding(emb1, colors=cd0[g, ], main=g, alpha=0.05)
}))


####### get one com, and recur
recur <- function(com) {
  com.sub <- lapply(levels(com), function(c) {
    print(c)
    cells <- names(com)[com==c]
    foo <- rownames(emb1) %in% cells
    names(foo) <- rownames(emb1)
    par(mfrow=c(1,1))
    plotEmbedding(emb1, groups=foo, show.legend=TRUE, main='PCA')

    cd <- normalizeCounts(pbmcA[, cells])
    class(cd)
    matnorm <- normalizeVariance(cd, plot=TRUE, details=TRUE, alpha=0.05)
    names(matnorm)

    nPcs = 30

    if(length(matnorm$ods)<=nPcs) {
      com <- rep(c, length(cells))
      names(com) <- cells
      return(list(com=com, model=NA))
    }
    mat <- log10(matnorm$mat+1)[matnorm$ods,]
    dim(mat)

    pcs <- fastPca(t(mat), nPcs=nPcs)
    m <- pcs$l
    colnames(m) <- paste0('PC', seq_len(nPcs))
    pcs.fin <- m

    emb <- Rtsne::Rtsne(as.matrix(pcs.fin), is_distance=FALSE, perplexity=10, verbose=TRUE, num_threads=4)$Y
    rownames(emb) <- rownames(pcs.fin)

    com <- getComMembership(t(as.matrix(pcs.fin)), k=30, method=igraph::cluster_walktrap)
    par(mfrow=c(1,2))
    plotEmbedding(emb, groups=com, show.legend=TRUE)
    plotEmbedding(emb1, groups=com, show.legend=TRUE)

    levels(com) <- paste0(c, '-', levels(com))

    ## each cell is its own group...can't be
    if(length(levels(com))==nrow(pcs.fin)) {
      com <- rep(c, length(cells))
      names(com) <- cells
      return(list(com=com, model=NA))
    }

    com.fin <- as.character(com)
    names(com.fin) <- names(com)

    if(length(unique(com.fin))==1) {
      return(list(com=com.fin, model=NA))
    } else {
      com.fin.stable <- getStableClusters(cm=t(cd), cols=factor(com.fin))
      if(length(unique(com.fin.stable))==1) {
        return(list(com=com.fin.stable, model=NA))
      }
      while(sum(as.character(com.fin)==as.character(com.fin.stable))!=length(com.fin.stable)) {
        com.fin <- com.fin.stable
        com.fin.stable <- getStableClusters(cm=t(cd), cols=factor(com.fin))
        if(length(unique(com.fin.stable))==1) {
          return(list(com=com.fin.stable, model=NA))
        }
      }

      plotEmbedding(emb, groups=com.fin.stable, show.legend=TRUE)
      plotEmbedding(emb1, groups=com.fin.stable, show.legend=TRUE)

      if(length(unique(com.fin.stable))==1) {
        return(list(com=com.fin.stable, model=NA))
      } else {
        model <- modelLda(mat, com.fin.stable)
        #reduction <- predict(model, data.frame(t(as.matrix(mat))))$x
        #ds <- getDifferentialGenes(cm=t(cd), cols=com.fin.stable)
        ##emb2 <- Rtsne::Rtsne(reduction, is_distance=FALSE, perplexity=30, verbose=TRUE, num_threads=4)$Y
        ##rownames(emb2) <- colnames(mat)
        #par(mfrow=c(1,3))
        #plotEmbedding(emb, groups=com.fin.stable, show.legend=TRUE, main="PCA")
        #plotEmbedding(emb2, groups=com.fin.stable, show.legend=TRUE, main="LDA")
        #plotEmbedding(emb1, groups=com.fin.stable, show.legend=TRUE)

        return(list(com=com.fin.stable, model=model))
      }
    }
  })
  return(com.sub)
}
com.sub <- recur(com0.stable)

com.fin <- unlist(lapply(com.sub, function(x) {
  foo <- as.character(x$com)
  names(foo) <- names(x$com)
  return(foo)
  }))
unique(com.fin)
com.fin <- factor(com.fin)

par(mfrow=c(1,3))
plotEmbedding(emb1, groups=com0.stable, show.legend=TRUE, main="Stable")
plotEmbedding(emb1, groups=com.fin, show.legend=TRUE, main="Recur")
plotEmbedding(emb1, groups=com.fin, show.legend=TRUE, shuffle.colors=TRUE)

models <- lapply(com.sub, function(x) x$model)
vi <- lapply(models, length) > 1
models <- models[vi]
models[[length(models)+1]] <- model # add back original
emb3 <- tsneLda(mat0, models, plot=FALSE)
par(mfrow=c(1,3))
plotEmbedding(emb3, groups=com0, show.legend=TRUE)
plotEmbedding(emb3, groups=com.fin, show.legend=TRUE)
plotEmbedding(emb3, groups=com.fin, shuffle.colors=TRUE)
