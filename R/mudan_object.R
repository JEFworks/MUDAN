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
        lds = NULL,
        com = NULL,
        model = NULL,
        emb = NULL,


        initialize =
            function(name = NA, counts = NA, verbose = TRUE, ncores=1, ...)
            {
                self$name <- name
                self$counts <- cleanCounts(counts, ...)
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
            function(reductionType='pcs', communityName="Infomap", communityMethod=igraph::cluster_infomap, k=30, minSize=5, ...)
            {
                com <- getComMembership(self[[reductionType]], method=communityMethod, k=k, ...)

                min.group.size <- minSize
                bad.groups <- names(which(table(com) < min.group.size))
                com[com %in% bad.groups] <- NA

                com <- factor(com)

                self$com[[reductionType]][[communityName]] <- com
            },

        modelCommunity =
            function(nGenes=min(1000, nrow(self$matnorm)*0.5), groups=NULL, communityName=NULL, ...)
            {
                vargenes <- getVariableGenes(self$matnorm, nGenes)
                ## test all coms
                if(is.null(groups)) {
                    if(is.null(communityName)) {
                        self$model <- lapply(self$com[['pcs']], function(c) {
                                                 model <- modelLda(mat=self$mat[vargenes,], com=c)
                                             })
                    } else {
                        self$model <- modelLda(self$mat[vargenes,], self$com[['pcs']][[communityName]], ...)
                    }
                } else {
                    self$model <- modelLda(self$mat[vargenes,], groups, ...)
                }
            },

        getMudanEmbedding =
            function(...)
            {
                results <- tsneLda(mat=self$mat, model=self$model, details=TRUE, ...)
                self$lds <- t(results$reduction)
                self$emb[['MUDAN']] <- results$emb
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
                    plotEmbedding(pcs.emb, ...)
                }

            },

        plot =
            function(reductionType, groups=NULL, communityName, embeddingType, ...)
            {
                if(is.null(groups)) {
                    plotEmbedding(self$emb[[embeddingType]], groups=self$com[[reductionType]][[communityName]], ...)
                } else {
                  plotEmbedding(self$emb[[embeddingType]], groups=groups, ...)
                }
            }
    )
)

