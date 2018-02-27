#' R6 object for storing all relevant output from MUDAN analysis in single entity
#'
#' @examples {
#' library(MUDAN)
#' data(pbmcA)
#' cd <- pbmcA
#' myMudanObject <- Mudan$new("PBMCA", cd)
#' myMudanObject$normalizeCounts()
#' myMudanObject$normalizeVariance(plot=TRUE)
#' myMudanObject$getPcs(nPcs=30, maxit=1000)
#' myMudanObject$getComMembership(communityName='Infomap',
#'      communityMethod=igraph::cluster_infomap, k=1/200*ncol(cd))
#' myMudanObject$getStableClusters(communityName='Infomap')
#' myMudanObject$modelCommunity(communityName='Infomap')
#' myMudanObject$getMudanEmbedding()
#' myMudanObject$getStandardEmbedding()
#' par(mfrow=c(1,2))
#' myMudanObject$plotEmbedding(communityName='Infomap', embeddingType='PCA',
#'      xlab=NULL, ylab=NULL, main='Standard')
#' myMudanObject$plotEmbedding(communityName='Infomap', embeddingType='MUDAN',
#'      xlab=NULL, ylab=NULL, main='MUDAN')
#' }
#'
#' @export
#'
Mudan <- R6::R6Class(
    "Mudan",
    public = list(
        name = NULL,
        cd = NULL,
        mat = NULL,
        matnorm = NULL,
        verbose = NULL,
        ncores = NULL,

        pcs = NULL,
        lds = NULL,
        com = NULL,
        model = NULL,
        preds = NULL,
        emb = NULL,

        gsf = NULL,
        ods = NULL,
        pv.sig = NULL,

        initialize =
            function(name = NA, counts = NA, verbose = TRUE, ncores=1, ...)
            {
                self$name <- name
                self$cd <- cleanCounts(counts, ...)
                self$verbose <- verbose
                self$ncores <- ncores
            },

        normalizeCounts =
            function(...)
            {
                self$mat <- normalizeCounts(self$cd, ...)
            },

        normalizeVariance =
            function(...)
            {
                matnorm.info <- normalizeVariance(self$mat, details=TRUE, ...) ## variance normalize
                ## remember normalization parameters
                gsf <- matnorm.info$df$gsf
                names(gsf) <- rownames(matnorm.info$df)
                self$gsf <- gsf
                self$ods <- rownames(matnorm.info$mat)[matnorm.info$ods]
                self$matnorm <- log10(matnorm.info$mat+1)
            },

        getPcs =
            function(nPcs=30, ...)
            {
                self$pcs <- getPcs(self$matnorm[self$ods,], nGenes=length(self$ods), nPcs=nPcs, ...)
            },

        getComMembership =
            function(communityName="Infomap", communityMethod=igraph::cluster_infomap, k=15, ...)
            {
                com <- getComMembership(self$pcs, method=communityMethod, k=k, ...)
                com <- factor(com)

                self$com[[communityName]][['ALL']] <- com
            },

        getStandardEmbedding =
          function(perplexity)
          {
            emb <- Rtsne::Rtsne(self$pcs, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y ## get tSNE embedding on PCs
            rownames(emb) <- rownames(self$pcs)
            self$emb[['PCA']] <- emb
          },

        getStableClusters =
          function(communityName="Infomap", min.group.size=10, min.diff.genes=10, z.threshold=3, ...)
          {
            ## get stable clusters
            comA <- getStableClusters(self$cd,
                                      self$com[[communityName]][['ALL']],
                                      self$matnorm,
                                      min.group.size=min.group.size,
                                      min.diff.genes=min.diff.genes,
                                      z.threshold=z.threshold,
                                      ...)
            self$com[[communityName]][['STABLE']] <- comA$com
            self$pv.sig[[communityName]] <- comA$pv.sig
          },

        modelCommunity =
            function(communityName="Infomap", groups=NULL, ...)
            {
              pv.sig.genes <- lapply(self$pv.sig[[communityName]], function(x) unique(unlist(lapply(x, rownames))))
              genes <- intersect(unique(unlist(pv.sig.genes)), rownames(matnorm))
              mn <- self$matnorm[genes,]

              if(is.null(groups)) {
                  self$model <- modelLda(mat=mn, com=self$com[[communityName]][['STABLE']], retest=FALSE, ...)
              } else {
                  self$model <- modelLda(mat=mn, com=groups, retest=FALSE, ...)
              }
            },

        getMudanEmbedding =
            function(perplexity=30, ...)
            {
                self$preds <- predict(self$model, data.frame(t(log10(as.matrix(self$mat[names(self$gsf),]*self$gsf+1)))))
                self$lds <- self$preds$x

                emb.lds <- Rtsne::Rtsne(self$lds, is_distance=FALSE, perplexity=perplexity, verbose=TRUE, ...)$Y ## get tSNE embedding on PCs
                rownames(emb.lds) <- rownames(self$lds)

                self$emb[['MUDAN']] <- emb.lds
            },

        plotEmbedding =
            function(groups=NULL, embeddingType, communityName, communityType='STABLE', ...)
            {
                if(is.null(groups)) {
                    plotEmbedding(self$emb[[embeddingType]], groups=self$com[[communityName]][[communityType]], ...)
                } else {
                    plotEmbedding(self$emb[[embeddingType]], groups=groups, ...)
                }
            }
    )
)

