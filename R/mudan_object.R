##' Mudan object
##'
##' @example
##' library(MUDAN)
##' source('../R/mudan.R')
##' data(pbmcA)
##' myMudanObject <- Mudan$new("PBMC", pbmcA)
##' myMudanObject$libSizeNormalize()
##' myMudanObject$varianceNormalize(plot=TRUE)
##' myMudanObject$dimensionalityReduction(maxit=1000)
##' myMudanObject$communityDetection(communityName='Infomap', communityMethod=igraph::cluster_infomap, k=10, minSize=10)
##' myMudanObject$communityDetection(communityName='Walktrap', communityMethod=igraph::cluster_walktrap, k=10, minSize=10)
##' myMudanObject$modelCommunity()
##' myMudanObject$getMudanEmbedding(plot=FALSE)
##' myMudanObject$getStandardEmbedding(plot=FALSE)
##' par(mfrow=c(2,2))
##' myMudanObject$plot(communityName='Infomap', embeddingType='MUDAN')
##' myMudanObject$plot(communityName='Walktrap', embeddingType='MUDAN')
##' myMudanObject$plot(communityName='Infomap', embeddingType='PCA')
##' myMudanObject$plot(communityName='Walktrap', embeddingType='PCA')
##'
Mudan <- R6::R6Class(
    "Mudan",
    public = list(
        name = NULL,
        counts = NULL,
        cd = NULL,
        mat = NULL,
        matnorm = NULL,
        verbose = NULL,
        ncores = NULL,

        pcs = NULL,
        com = NULL,
        model = NULL,
        emb = NULL,

        initialize =
            function(name = NA, counts = NA, verbose = TRUE, ncores=10, ...)
            {
                self$name <- name
                self$counts <- cleanCounts(as.matrix(counts), ...)
                self$verbose <- verbose
                self$ncores <- ncores
            },

        libSizeNormalize =
            function(...)
            {
                self$cd <- normalizeCounts(self$counts, ...)
                self$mat <- log10(self$cd+1)
            },

        varianceNormalize =
            function(...)
            {
                self$matnorm <- log10(normalizeVariance(self$cd, ...)+1)
            },

        dimensionalityReduction =
            function(...)
            {
                self$pcs <- getPcs(self$matnorm, ...)
            },

        communityDetection =
            function(communityName, communityMethod, minSize=5, ...)
            {
                com <- getComMembership(self$pcs, method=communityMethod, ...)

                min.group.size <- minSize
                bad.groups <- names(which(table(com) < min.group.size))
                com[com %in% bad.groups] <- NA

                self$com[[communityName]] <- com
            },

        modelCommunity =
            function(nGenes=1000, communityName=1, ...)
            {
                vargenes <- getVariableGenes(self$matnorm, nGenes)
                self$model <- modelLda(self$mat[vargenes,], self$com[[communityName]], ...)
            },

        getMudanEmbedding =
            function(communityName = 1, ...)
            {
                self$emb[['MUDAN']] <- tsneLda(mat=self$mat, model=self$model, com=self$com[[communityName]], ...)
            },

        getStandardEmbedding =
            function(plot=TRUE, do.par=TRUE, ...)
            {
                pcs.emb <- Rtsne::Rtsne(t(self$pcs), is_distance=FALSE, num_threads=self$ncores)$Y
                rownames(pcs.emb) <- colnames(self$pcs)
                self$emb[['PCA']] <- pcs.emb
                if(plot) {
                    if(do.par) {
                        par(mar = c(0.5,0.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
                    }
                    plotEmbedding(pcs.emb, groups=self$com, ...)
                }

            },

        plot =
            function(communityName, embeddingType, ...)
            {
                plotEmbedding(self$emb[[embeddingType]], groups=self$com[[communityName]], ...)
            }
    )
)

