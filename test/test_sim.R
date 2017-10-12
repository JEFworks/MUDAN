library(mudan)
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
test.sim <- function() {
    ## simulate data
    data <- simulate.data()
    mat <- data$mat
    group <- data$group

    ## discover subpopulations, train discriminating model
    cd <- as.matrix(mat)
    dim(cd)
    matnorm <- normalizeVariance(mat=cd, plot=TRUE, details=FALSE)
    pcs <- getPcs(matnorm, nGenes=100, nPcs=10)
    models <- lapply(c(15, 30, 50, 100, 150), function(k) {
        com <- getKnnMembership(pcs, k=k, method=igraph::cluster_infomap)
        reference.model <- modelLda(matnorm, com)
    })

    ## apply model to any normalization common to other datasets
    mat <- counts2cpms(cd)
    reference.emb <- tsneLda(mat=mat, model=models, com=group, perplexity=30, plot=FALSE)

    ## plot resulting embedding
    plotEmbedding(reference.emb,
                  groups=group,
                  mark.clusters = FALSE,
                  show.legend=TRUE,
                  legend.x = 'bottomleft',
                  main="MUDAN with true cell annotations")
}

