##################### Testing MUDAN using simulations
library(MUDAN)
source('../R/mudan.R')

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
data <- simulate.data(G=3, M=5000, N=1000)
mat <- cleanCounts(data$mat)
group <- data$group

## discover subpopulations, train discriminating model
cd <- normalizeCounts(mat)
matnorm <- log10(normalizeVariance(cd, plot=TRUE, details=FALSE)+1)
pcs <- getPcs(matnorm, nGenes=100, nPcs=10)
com <- getKnnMembership(pcs, k=30, method=igraph::cluster_infomap)
model <- modelLda(matnorm, com)
emb <- tsneLda(mat=mat, model=model, com=com, perplexity=30, details=FALSE, plot=FALSE)

## diff genes
diffGenes <- getDifferentialGenes(mat, com)
## markers
markers <- unlist(lapply(diffGenes$info, function(x) {
    rownames(x)[x$marker.auc>0.9]
}))

## plot resulting embedding
heatmap(mat[, names(sort(com))], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(length(levels(group)))[group[names(sort(com))]], labCol=FALSE, labRow=FALSE, main='\nSimulated Data')
par(mfrow=c(2,2), mar=rep(5,4))
plotEmbedding(emb,
              groups=group,
              mark.clusters = FALSE,
              show.legend=TRUE,
              legend.x = 'bottomleft',
              main="MUDAN with true\nsimulated cell annotations")
plotEmbedding(emb,
              groups=com,
              mark.clusters = FALSE,
              show.legend=TRUE,
              legend.x = 'bottomleft',
              main="MUDAN with detected\ninfomap community annotations")
heatmap(mat[markers,names(sort(com))], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", labCol=FALSE, labRow=FALSE, main='\nIdentified markers: AUC > 0.9', ColSideColors=rainbow(length(levels(com)))[com[names(sort(com))]])


