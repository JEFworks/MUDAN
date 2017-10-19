############################### Using MUDAN for projections and deconvolutions
library(MUDAN)
source('../R/mudan.R')

###############################
## Projection
###############################

############################### train on annotated pbmcs
library(Matrix)
data("referenceAnnot")
data("referenceCounts")

cd <- normalizeCounts(referenceCounts)
reference.mat <- mat <- log10(cd+1)
matnorm <- log10(normalizeVariance(cd, plot=TRUE, details=FALSE)+1)
pcs <- getPcs(matnorm, maxit=1000)
pcs.emb <- Rtsne::Rtsne(t(pcs), is_distance=FALSE, perplexity=100, num_threads=10)$Y
rownames(pcs.emb) <- colnames(matnorm)

par(mfrow=c(1,1), mar=rep(5,4))
plotEmbedding(pcs.emb,
              groups=referenceAnnot,
              mark.clusters = TRUE,
              show.legend=TRUE,
              legend.x = 'bottomleft',
              main="Reference")

############################### get markers
library(parallel)
## significant differential expression only
diffGenes <- getDifferentialGenes(mat, referenceAnnot)
## also useful in distinguishing
markers <- lapply(diffGenes$info, function(x) {
    x <- x[x$percent.expressing>0.5,]
    x <- x[x$log2.fold.change>0.5,]
    x <- x[x$marker.auc>0.5,]
    print(x)
    rownames(x)
})
## check markers per group; make sure > 0; more seems to be better
lapply(markers, length)
markers <- unique(unlist(markers))
length(markers)

############################## predict on new unannotated pbmc samples
## helper function
testPbmcs <- function(pbmcA) {
    cd <- normalizeCounts(cleanCounts(pbmcA))
    mat <- log10(cd+1)
    matnorm <- log10(normalizeVariance(cd, plot=TRUE, details=FALSE)+1)

    ## make prediction model
    markers <- intersect(markers, rownames(mat))
    print(length(markers))
    pred <- modelLda(reference.mat[markers,], referenceAnnot)

    ## predict
    r1 <- predict(pred, data.frame(t(as.matrix(mat))))
    r1pred <- r1$class
    names(r1pred) <- colnames(mat)

    ## plot
    pcs <- getPcs(matnorm, maxit=1000)
    pcs.emb <- Rtsne::Rtsne(t(pcs), is_distance=FALSE, perplexity=100, num_threads=10)$Y
    rownames(pcs.emb) <- colnames(matnorm)
    par(mfrow=c(1,1))
    plotEmbedding(pcs.emb,
                  groups=r1pred,
                  mark.clusters = TRUE,
                  show.legend=TRUE,
                  legend.x = 'bottomleft',
                  alpha=0.1)

    return(list(pcs=pcs, pcs.emb=pcs.emb, r1pred=r1pred))
}

## patient A
data(pbmcA)
pbmcA.results <- testPbmcs(as.matrix(pbmcA))
## patient B
data(pbmcB)
pbmcB.results <- testPbmcs(as.matrix(pbmcB))
## patient C
data(pbmcC)
pbmcC.results <- testPbmcs(as.matrix(pbmcC))

refine <- function(pbmcA.results) {
    ## refine based on knn
    mat <- pbmcA.results$pcs
    com <- pbmcA.results$r1pred

    ## model based refinement until convergence
    pf <- 1
    i <- 1
    while(pf > 0.01) {
        print(paste0('Refine iteration: ', i))
        df <- data.frame(celltype=com, t(mat))
        model <- MASS::lda(celltype ~ ., data=df)
        model.output <- predict(model, df)
        com.refine <- model.output$class
        names(com.refine) <- rownames(df)

        pf <- table(com.refine==com)[[1]]/length(com) ## percent false
        print(pf)

        com <- com.refine
    }

    par(mfrow=c(1,2))
    plotEmbedding(pbmcA.results$pcs.emb,
                  groups=pbmcA.results$r1pred,
                  mark.clusters = TRUE,
                  show.legend=TRUE,
                  legend.x = 'topleft',
                  alpha=0.1, main="Patient A with Predicted Annotations", cex.main=0.8)
    plotEmbedding(pbmcA.results$pcs.emb,
                  groups=com.refine,
                  mark.clusters = TRUE,
                  show.legend=TRUE,
                  legend.x = 'topleft',
                  alpha=0.1, main="Patient A with Refined Predicted Annotations", cex.main=0.8)

    return(com.refine)
}

pbmcA.refine <- refine(pbmcA.results)
pbmcB.refine <- refine(pbmcB.results)
pbmcC.refine <- refine(pbmcC.results)

png('../docs/tester_project.png', width=1200, height=300)
par(mfrow=c(1,4))
plotEmbedding(pcs.emb,
              groups=referenceAnnot,
              mark.clusters = TRUE,
              show.legend=TRUE,
              legend.x = 'topleft',
              main="Reference with True Annotations", cex.main=0.8)
plotEmbedding(pbmcA.results$pcs.emb,
              groups=pbmcA.refine,
              mark.clusters = TRUE,
              show.legend=TRUE,
              legend.x = 'topleft',
              alpha=0.1, main="Patient A with Predicted Annotations", cex.main=0.8)
plotEmbedding(pbmcB.results$pcs.emb,
              groups=pbmcB.refine,
              mark.clusters = TRUE,
              show.legend=TRUE,
              legend.x = 'topleft',
              alpha=0.1, main="Patient B with Predicted Annotations", cex.main=0.8)
plotEmbedding(pbmcC.results$pcs.emb,
              groups=pbmcC.refine,
              mark.clusters = TRUE,
              show.legend=TRUE,
              legend.x = 'topleft',
              alpha=0.1, main="Patient C with Predicted Annotations", cex.main=0.8)
dev.off()

## look at proportion predicted
table(pbmcA.refine)/length(pbmcA.refine)
table(pbmcB.refine)/length(pbmcB.refine)
table(pbmcC.refine)/length(pbmcC.refine)


###############################
## Deconvolve
###############################

## make artificial bulks
reference.bulk <- rowSums(referenceCounts)
pbmcA.bulk <- rowSums(pbmcA)
pbmcB.bulk <- rowSums(pbmcB)
pbmcC.bulk <- rowSums(pbmcC)

train <- do.call(cbind, lapply(levels(factor(referenceAnnot)), function(g) {
    cells <- names(referenceAnnot)[referenceAnnot==g]
    rowSums(referenceCounts[, cells])
}))
colnames(train) <- levels(factor(referenceAnnot))
train <- t(t(train)/colSums(train))*1e6
train <- train[markers,]
head(train)

test <- list(
    reference=(reference.bulk)/sum(reference.bulk)*1e6,
    pbmcA=(pbmcA.bulk)/sum(pbmcA.bulk)*1e6,
    pbmcB=(pbmcB.bulk)/sum(pbmcB.bulk)*1e6,
    pbmcC=(pbmcC.bulk)/sum(pbmcC.bulk)*1e6
)

## assess proportions in mixtures
library(nnls) ## non negative least square regression
results <-
    lapply(test, function(y) {
        fit <- nnls(A=train[markers,], b=y[markers])

        ## permutation to assess significance
        rand <- do.call(cbind, lapply(1:1000, function(seed) {
                                   set.seed(seed)
                                   y.reorder <- sample(y)
                                   names(y.reorder) <- names(y)
                                   fit <- nnls(A=train[markers,], b=y.reorder[markers])
                                   fit$x
                               }))
        z.score <- (fit$x-apply(rand, 1, mean))/apply(rand, 1, sd)
        z.score[is.nan(z.score)] <- 0

        result <- list(prop=fit$x, z=z.score)

        return(result)
    })
results.prop <- do.call(cbind, lapply(results, function(x) x$prop))
colnames(results.prop) <- names(test)
rownames(results.prop) <- colnames(train)
results.prop

results.z <- do.call(rbind, lapply(results, function(x) x$z))
colnames(results.z) <- colnames(train)
results.z

## resulting proportion of cells in mixture
proportion <- t(results.prop)/colSums(results.prop)
proportion

## truth
truth <- rbind(
    reference=table(referenceAnnot)/length(referenceAnnot),
    pbmcA=table(pbmcA.refine)/length(pbmcA.refine),
    pbmcB=table(pbmcB.refine)/length(pbmcB.refine),
    pbmcC=table(pbmcC.refine)/length(pbmcC.refine)
)
truth

results.z.col <- results.z
results.z.col[results.z.col<0] <- 0
results.z.col <- round(results.z.col*99)+1
results.z.col[results.z.col>100] <- 100

par(mfrow=c(1,4), mar=rep(5,4))
plot(truth[1,], proportion[1,], col=colorRampPalette(c('black', 'red'))(100)[results.z.col[1,]], pch=16, cex=2)
plot(truth[2,], proportion[2,], col=colorRampPalette(c('black', 'red'))(100)[results.z.col[2,]], pch=16, cex=2)
plot(truth[3,], proportion[3,], col=colorRampPalette(c('black', 'red'))(100)[results.z.col[3,]], pch=16, cex=2)
plot(truth[4,], proportion[4,], col=colorRampPalette(c('black', 'red'))(100)[results.z.col[4,]], pch=16, cex=2)
