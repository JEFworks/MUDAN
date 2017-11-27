##' Filter counts matrix
##'
##' Filter counts matrix based on gene and cell requirements
##'
##' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
##' @param min.lib.size Minimum number of genes detected in a cell. Cells with fewer genes will be removed (default: 1000)
##' @param max.lib.size Maximum number of genes detected in a cell. Cells with more genes will be removed (default: 8000)
##' @param min.reads Minimum number of reads per gene. Genes with fewer reads will be removed (default: 10)
##' @param min.detected Minimum number of cells a gene must be seen in. Genes not seen in a sufficient number of cells will be removed (default: 5)
##' @param verbose Verbosity (default: TRUE)
##'
##' @return a filtered read count matrix
##'
##' @examples {
##' data(pbmcA)
##' dim(pbmcA)
##' mat <- cleanCounts(pbmcA)
##' dim(mat)
##' }
##'
##' @export
##'
cleanCounts <- function(counts, min.lib.size = 1000, max.lib.size = 8000, min.reads = 1, min.detected = 1, verbose=TRUE) {
    if(class(counts)!='matrix') {
        counts <- as.matrix(counts)
    }

    if(verbose) {
        print(paste0('Filtering matrix with ', ncol(counts), ' cells and ', nrow(counts), ' genes ...'))
    }

    ## filter out low-gene cells
    counts <- counts[, Rfast::colsums(counts)>min.lib.size]
    ## filter out potential doublets
    counts <- counts[, Rfast::colsums(counts)<max.lib.size]
    ## remove genes that don't have many reads
    counts <- counts[Rfast::rowsums(counts)>min.reads, ]
    ## remove genes that are not seen in a sufficient number of cells
    counts <- counts[Rfast::rowsums(counts>0)>min.detected, ]

    if(verbose) {
        print(paste0('Resulting matrix has ', ncol(counts), ' cells and ', nrow(counts), ' genes'))
    }

    return(counts)
}


##' Normalizes counts to CPM
##'
##' Normalizes raw counts to log10 counts per million with pseudocount
##'
##' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
##' @param pseudocount Pseudocount for log transformation. (default: 1)
##' @param verbose Verbosity (default: TRUE)
##'
##' @return a normalized matrix
##'
##' @examples {
##' data(pbmcA)
##' mat <- counts2cpms(pbmcA)
##' }
##'
##' @export
##'
normalizeCounts <- function(counts, depthScale=1e6, verbose=TRUE) {
    if(class(counts)!='matrix') {
        counts <- as.matrix(counts)
    }
    if(verbose) {
        print(paste0('Normalizing matrix with ', ncol(counts), ' cells and ', nrow(counts), ' genes'))
    }

    counts <- t(t(counts)/Rfast::colsums(counts))
    counts <- counts*depthScale

    return(counts)
}


##' Normalize gene expression variance relative to transcriptome-wide expectations
##' (Modified from SCDE/PAGODA2 code; now uses Rfast)
##'
##' Normalizes gene expression magnitudes to with respect to its ratio to the
##' transcriptome-wide expectation as determined by local regression on expression magnitude
##'
##' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
##' @param gam.k Generalized additive model parameter; the dimension of the basis used to represent the smooth term (default: 5)
##' @param alpha Significance threshold (default: 0.05)
##' @param plot Whether to plot the results (default: FALSE)
##' @param use.unadjusted.pvals If true, will apply BH correction (default: FALSE)
##' @param do.par Whether to adjust par for plotting if plotting (default: TRUE)
##' @param max.adjusted.variance Ceiling on maximum variance after normalization to prevent infinites (default: 1e3)
##' @param min.adjusted.variance Floor on minimum variance after normalization (default: 1e-3)
##' @param verbose Verbosity (default: TRUE)
##' @param details If true, will return data frame of normalization parameters. Else will return normalized matrix.(default: FALSE)
##'
##' @return If details is true, will return data frame of normalization parameters. Else will return normalized matrix.
##'
##' @examples {
##' data(pbmcA)
##' mat <- cleanCounts(pbmcA)
##' mat <- normalizeVariance(mat)
##' }
##'
##' @importFrom mgcv s
##'
##' @export
##'
normalizeVariance <- function(cd, gam.k=5, alpha=0.05, plot=FALSE, use.unadjusted.pvals=FALSE, do.par=TRUE, max.adjusted.variance=1e3, min.adjusted.variance=1e-3, verbose=TRUE, details=FALSE) {
    if(class(cd)!='matrix') {
        cd <- as.matrix(cd)
    }

    mat <- t(cd) ## make rows as cells, cols as genes

    if(verbose) {
        print("Calculating variance fit ...")
    }
    dfm <- log(Rfast::colmeans(mat))
    dfv <- log(Rfast::colVars(mat))
    names(dfm) <- names(dfv) <- colnames(mat)
    df <- data.frame(m=dfm, v=dfv)

    vi <- which(is.finite(dfv))

    if(length(vi)<gam.k*1.5) { gam.k=1 } ## too few genes

    if(gam.k<2) {
        if(verbose) {
            print("Using lm ...")
        }
        m <- lm(v ~ m, data = df[vi,])
    } else {
        if(verbose) {
            print(paste0("Using gam with k=", gam.k, "..."))
        }
        fm <- as.formula(sprintf("v ~ s(m, k = %s)", gam.k))
        m <- mgcv::gam(fm, data = df[vi,])
    }
    df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
    n.cells <- ncol(mat)
    n.obs <- nrow(mat)
    df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=T))
    df$lpa <- bh.adjust(df$lp,log=TRUE)
    df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)

    if(use.unadjusted.pvals) {
        ods <- which(df$lp<log(alpha))
    } else {
        ods <- which(df$lpa<log(alpha))
    }
    if(verbose) {
        print(paste0(length(ods), ' overdispersed genes ... ' ))
    }

    df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
    df$gsf[!is.finite(df$gsf)] <- 0;

    if(plot) {
        if(do.par) {
            par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
        }
        smoothScatter(df$m,df$v,main='',xlab='log10[ magnitude ]',ylab='log10[ variance ]')
        grid <- seq(min(df$m[vi]),max(df$m[vi]),length.out=1000)
        lines(grid,predict(m,newdata=data.frame(m=grid)),col="blue")
        if(length(ods)>0) {
            points(df$m[ods],df$v[ods],pch='.',col=2,cex=1)
        }
        smoothScatter(df$m[vi],df$qv[vi],xlab='log10[ magnitude ]',ylab='',main='adjusted')
        abline(h=1,lty=2,col=8)
        if(is.finite(max.adjusted.variance)) { abline(h=max.adjusted.variance,lty=2,col=1) }
        points(df$m[ods],df$qv[ods],col=2,pch='.')
    }

    if(!details) {
        ## variance normalize
        norm.mat <- cd*df$gsf
        return(norm.mat)
    }
    else {
        ## return normalization factor
        return(df)
    }
}
## BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
    nai <- which(!is.na(x))
    ox <- x
    x <- x[nai]
    id <- order(x, decreasing = FALSE)
    if(log) {
        q <- x[id] + log(length(x)/seq_along(x))
    } else {
        q <- x[id]*length(x)/seq_along(x)
    }
    a <- rev(cummin(rev(q)))[order(id)]
    ox[nai] <- a
    ox
}


##' Dimensionality reduction by PCA
##'
##' Dimensionality reduction using PCA by computing principal components using highly variable genes
##'
##' @param mat Variance normalized gene expression matrix.
##' @param nGenes Number of most variable genes. (default: 1000)
##' @param nPcs Number of principal components. (default: 100)
##' @param verbose Verbosity (default: TRUE)
##'
##' @return Matrix with columns as cells and rows as principal component eigenvectors.
##'
##' @examples {
##' data(pbmcA)
##' mat <- cleanCounts(pbmcA)
##' mat <- normalizeVariance(mat)
##' pcs <- getPcs(mat)
##' }
##'
##' @export
##'
getPcs <- function(mat, nGenes = min(nrow(mat), 1000), nPcs = 100, verbose=TRUE, ...) {
    if(class(mat)!='matrix') {
        mat <- as.matrix(mat)
    }

    if(verbose) {
        print(paste0('Identifying top ', nPcs, ' PCs on ', nGenes, ' most variable genes ...'))
    }

    ## get variable genes
    vgenes <- getVariableGenes(mat, nGenes)
    mat <- mat[vgenes,]

    ## get PCs
    pcs <- fastPca(t(mat), nPcs, ...)
    m <- t(pcs$l)
    colnames(m) <- colnames(mat)
    rownames(m) <- paste0('PC', seq_len(nPcs))

    return(m)
}
fastPca <- function(m, nPcs=2, tol=1e-10, scale=FALSE, center=FALSE, transpose=FALSE, ...) {
    ## note transpose is meant to speed up calculations when neither scaling nor centering is required
    if(transpose) {
        if(center) {
            m <- m-Rfast::rowmeans(m)
        }
        if(scale) {
            m <- m/sqrt(Rfast::rowsums(m*m))
        }
        a <- irlba::irlba(tcrossprod(m)/(ncol(m)-1), nu=0, nv=nPcs, tol=tol, ...)
        a$l <- t(t(a$v) %*% m)
    } else {
        if(scale||center) {
            m <- scale(m, scale=scale, center=center)
        }
        a <- irlba::irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs, tol=tol, ...)
        a$l <- m %*% a$v
    }
    return(a)
}
getVariableGenes <- function(mat, nGenes) {
    vi <- Rfast::rowVars(mat)
    names(vi) <- rownames(mat)
    vgenes <- names(sort(vi, decreasing=TRUE)[seq_len(nGenes)])
    return(vgenes)
}

##' Get clusters by community detection on approximate nearest neighbors
##'
##' Group cells into clusters based on graph-based community detection on approximate nearest neighbors
##'
##' @param mat Matrix of cells as columns. Features as rows (such as PCs).
##' @param k K-nearest neighbor parameter.
##' @param method Community detection method from igraph. (default: igraph::cluster_walktrap)
##' @param verbose Verbosity (default: TRUE)
##'
##' @return Vector of community annotations
##'
##' @examples {
##' data(pbmcA)
##' mat <- cleanCounts(pbmcA)
##' mat <- normalizeVariance(mat)
##' pcs <- getPcs(mat)
##' com <- getComMembership(pcs, k=30)
##' }
##'
##' @export
##'
getComMembership <- function(mat, k, method=igraph::cluster_walktrap, verbose=TRUE, details=FALSE) {
    if(verbose) {
        print("finding approximate nearest neighbors ...")
    }
    knn <- RANN::nn2(t(mat), k=k)[[1]]

    ## convert to adjacency matrix
    adj <- matrix(0, ncol(mat), ncol(mat))
    rownames(adj) <- colnames(adj) <- colnames(mat)
    invisible(lapply(seq_len(ncol(mat)), function(i) {
        adj[i,colnames(mat)[knn[i,]]] <<- 1
    }))

    ## convert to graph for clustering
    if(verbose) {
        print("calculating clustering ...")
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- igraph::simplify(g)
    km <- method(g)
    if(verbose) {
        mod <- igraph::modularity(km)
        if(mod < 0.3) { print('WARNING') }
        print(paste0('graph modularity: ', mod))
    }

    ## community membership
    com <- km$membership
    names(com) <- km$names
    com <- factor(com)

    if(verbose) {
        print("identifying cluster membership ...")
        print(table(com))
    }
    if(details) {
        return(list(com=com, mod=mod, g=g))
    }
    else {
        return(com)
    }
}
## Test a set of ks for modularity
## TODO: Modularity is a function of k; need permutation to establish baseline and see if observed modularity for a particular k is significantly better than random
optimizeModularity <- function(mat, ks=c(15,30,50,100), method=igraph::cluster_walktrap, verbose=TRUE) {
    results <- lapply(ks, function(k) {
        if(verbose) {
            print(paste0('testing k:', k, ' ...'))
        }
        getComMembership(mat, k, details=TRUE)
    })
    mods <- sapply(results, function(x) x$mod)
    i <- which(mods==max(mods))
    if(verbose) {
        print(paste0('optimal k:', ks[i], ' ...'))
    }
    com.final <- results[[i]]$com
    return(com.final)
}
## When there are too many cells, getKNNMembership runs into memory issues
## Also just takes too long
## Need subsampling
getApproxComMembership <- function(mat, k, nsubsample=ncol(mat)*0.5, method=igraph::cluster_walktrap, seed=0, verbose=TRUE) {

    if(verbose) {
        print(paste0('Subsampling from ', ncol(mat), ' cells to ', nsubsample, ' ... '))
    }
    ## random subsampling
    ## TODO: density based downsampling
    set.seed(seed)
    subsample <- sample(colnames(mat), nsubsample)

    if(verbose) {
        print('Identifying cluster membership for subsample ... ')
    }
    pcs.sub <- mat[, subsample]
    com.sub <- getComMembership(pcs.sub, k=k, method=method)

    if(verbose) {
        print('Imputing cluster membership for rest of cells ... ')
    }

    ## Use neighbor voting
    ## data <- mat[, subsample]
    ## query <- mat[, setdiff(colnames(mat), subsample)]
    ## knn <- RANN::nn2(t(data), t(query), k=k2)[[1]]
    ## rownames(knn) <- colnames(query)
    ## com.nonsub <- unlist(apply(knn, 1, function(x) {
    ##         ## nearest neighbors in data
    ##         nn <- colnames(data)[x]
    ##         ## look at their cell type annotations
    ##         nn.com <- com.sub[nn]
    ##         ## get most frequent annotation
    ##         return(names(sort(table(nn.com), decreasing=TRUE)[1]))
    ## }))
    ## com.all <- factor(c(com.sub, com.nonsub)[colnames(mat)])

    ## Use model instead
    ## Inspired by DenSVM
    df.sub <- data.frame(celltype=com.sub, t(pcs.sub))
    model <- MASS::lda(celltype ~ ., data=df.sub)
    df.all <- data.frame(t(mat))
    model.output <- predict(model, df.all)
    com.all <- model.output$class
    names(com.all) <- rownames(df.all)

    if(verbose) {
        print("LDA prediction accuracy for subsample ...")
        print(table(com.all[names(com.sub)]==com.sub))
    }

    return(com.all)
}


##' Linear discriminant analysis model
##'
##' Identifies components that maximally discriminate among groups using a linear discriminant analysis model
##'
##' @param mat Expression matrix with cells as columns, transferable features such as genes as rows.
##' @param com Community annotations
##' @param ncores Number of cores for parallelization.
##' @param verbose Verbosity (default: TRUE)
##' @param retest Whether to retest model for accuracy
##'
##' @return LDA model
##'
##' @examples {
##' data(pbmcA)
##' mat <- cleanCounts(pbmcA)
##' mat <- normalizeVariance(mat)
##' pcs <- getPcs(mat)
##' com <- getKnnMembership(pcs, k=30)
##' model <- modelLda(mat, com)
##' }
##'
##' @export
##'
modelLda.multicore <- function(mat, com, ncores=10, verbose=TRUE, retest=TRUE) {
    ## split features randomly into ncores groups
    ## ISSUE: each set of LDs only dependent on genes in group...not recommended
    set.seed(0)
    feature.set <- split(sample(rownames(mat)), seq_len(ncores), drop = FALSE)

    ## train model for each feature set in parallele to speed things up
    models <- parallel::mclapply(feature.set, function(vi) {
        ## make data frame
        df <- data.frame(celltype=com, as.matrix(t(mat[vi,])))

        if(verbose) {
            print("calculating LDA ...")
        }
        model <- MASS::lda(celltype ~ ., data=df)

        if(retest) {
            ## predict our data based on model
            model.output <- predict(model, df)
            if(verbose) {
                print("LDA prediction accuracy ...")
                print(table(model.output$class==com))
            }
        }

        return(model)
    }, mc.cores=ncores)

    return(models)
}
modelLda <- function(mat, com, nfeatures=nrow(mat), random=FALSE, verbose=TRUE, retest=TRUE) {
    ## filter to reduce feature modeling space
    if(nfeatures < nrow(mat)) {
        if(random) {
            ## random
            set.seed(0)
            vi <- sample(rownames(mat), nfeatures)
            mat <- mat[vi,]
        } else {
            ## most variable (assumes properly variance normalized matrix...fix)
            vi <- getVariableGenes(mat, nfeatures)
            mat <- mat[vi,]
        }
    }

    df <- data.frame(celltype=com, as.matrix(t(mat)))

    if(verbose) {
        print("calculating LDA ...")
    }
    model <- MASS::lda(celltype ~ ., data=df)

    if(retest) {
        ## predict our data based on model
        model.output <- predict(model, df)
        if(verbose) {
            print("LDA prediction accuracy ...")
            print(table(model.output$class==com))
        }
    }

    return(model)
}



getClusterGeneInfo <- function(mat, com, verbose=FALSE, ncores=10) {

    if(class(mat)!='matrix') {
        mat <- as.matrix(mat)
    }
    if(class(com)!='factor') {
        com <- factor(com)
    }
    if(!all(colnames(mat) %in% names(com))) { warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")}
    ## determine a subset of cells that's in the cols and cols[cell]!=NA
    valid.cells <- colnames(mat) %in% names(com)[!is.na(com)]
    if(!all(valid.cells)) {
        ## take a subset of the count matrix
        mat <- mat[, valid.cells]
    }
    ## reorder cols
    cols <- com
    cols <- as.factor(cols[match(colnames(mat),names(cols))])
    cols <- as.factor(cols)

    ## run simple wilcoxon test comparing each group with the rest
    diffgenes <- lapply(levels(cols), function(g) {
        if(verbose) {
            print(paste0("Testing group ", g, " ... "))
        }

        x <- mat[,cols==g]
        y <- mat[,cols!=g]

        ## aucs and pvalues
        aucs.pvs <- parallel::mclapply(seq_len(nrow(x)), function(i) {
            auc <- markerAuc(x[i,], y[i,])
            pv <- wilcox.test(x[i,], y[i,], alternative="greater")$p.value
            return(list(auc, pv))
        }, mc.cores=ncores)
        aucs <- sapply(aucs.pvs, function(x) x[[1]])
        pvs <- sapply(aucs.pvs, function(x) x[[2]])

        ## assess average fold change
        xm <- Rfast::rowmeans(x)
        ym <- Rfast::rowmeans(y)
        fc <- log2(xm/ym)
        ## percent expressing
        pe <- Rfast::rowsums(x>0)/ncol(x)
        ## z-score
        zs <- abs(qnorm(1 - (pvs/2), lower.tail=FALSE))*sign(fc)

        ## store p value
        df <- data.frame(p.value=pvs, marker.auc=aucs, log2.fold.change=fc, percent.expressing=pe, z.score=zs)
        rownames(df) <- rownames(mat)

        return(df)
    })
    names(diffgenes) <- levels(cols)

    return(diffgenes)
}
## Test how well does a gene discriminates a population
## adapted from seurat package
markerAuc <- function(x, y) {
    pred <- ROCR::prediction(
        predictions = c(x, y),
        labels = c(rep(x = 1, length(x)), rep(x = 0, length(y))),
        label.ordering = 0:1
    )
    perf <- ROCR::performance(pred, measure = "auc")
    auc <- perf@y.values[[1]]
    return(auc)
}

getDifferentialGenes <- function(cm, cols, verbose=TRUE) {

  ## match matrix rownames (cells) and group annotations
  if(!all(rownames(cm) %in% names(cols))) { warning("Cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")}
  ## determine a subset of cells that's in the cols and cols[cell]!=NA
  valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)];
  if(!all(valid.cells)) {
    ## take a subset of the count matrix
    cm <- cm[valid.cells,]
  }
  ## reorder cols
  cols <- as.factor(cols[match(rownames(cm),names(cols))]);
  cols <- as.factor(cols);

  if(verbose) {
    print(paste0("Running differential expression with ",length(levels(cols))," clusters ... "))
  }

  ## modified from pagoda2
  # run wilcoxon test comparing each group with the rest
  # calculate rank per column (per-gene) average rank matrix
  xr <- apply(cm, 2, function(foo) {
    #foo[foo==0] <- NA
    bar <- Rfast::Rank(foo)
    #bar[is.na(foo)] <- 0
    bar[foo==0] <- 0
    bar
  }); rownames(xr) <- rownames(cm)
  range(xr[,1])
  xr[1:5,1:5]
  cm[1:5,1:5]
  ##xr <- sparse_matrix_column_ranks(cm);

  # calculate rank sums per group
  grs <- do.call(rbind, lapply(levels(cols), function(g) Rfast::colsums(xr[cols==g,])))
  rownames(grs) <- levels(cols); colnames(grs) <- colnames(xr)
  grs[1:5,1:5]
  ##grs <- colSumByFac(xr,as.integer(cols))[-1,,drop=F]

  # calculate number of non-zero entries per group
  gnzz <- do.call(rbind, lapply(levels(cols), function(g) Rfast::colsums(xr[cols==g,]>0)))
  rownames(gnzz) <- levels(cols); colnames(gnzz) <- colnames(xr)
  gnzz[1:5,1:5]
  #xr@x <- numeric(length(xr@x))+1
  #gnzz <- colSumByFac(xr,as.integer(cols))[-1,,drop=F]

  #group.size <- as.numeric(tapply(cols,cols,length));
  #group.size <- as.numeric(tapply(cols,cols,length))[1:nrow(gnzz)]; group.size[is.na(group.size)]<-0; # trailing empty levels are cut off by colSumByFac
  group.size <- as.numeric(table(cols))

  # add contribution of zero entries to the grs
  gnz <- (group.size-gnzz)

  # rank of a 0 entry for each gene
  #zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
  zero.ranks <- apply(cm, 2, function(foo) {
    bar <- Rfast::Rank(foo)
    bar[foo==0][1]
  })
  ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2

  # standardize
  n1n2 <- group.size*(nrow(cm)-group.size);
  # usigma <- sqrt(n1n2*(nrow(cm)+1)/12) # without tie correction
  # correcting for 0 ties, of which there are plenty
  #usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
  usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
  # standardized U value- z score
  x <- t((ustat - n1n2/2)/usigma);

  # correct for multiple hypothesis
  x <- matrix(qnorm(bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE,
                                    log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),
              ncol = ncol(x)) * sign(x)
  rownames(x) <- colnames(cm)
  colnames(x) <- levels(cols)[1:ncol(x)]

  # add fold change information
  log.gene.av <- log2(Rfast::colmeans(cm));
  group.gene.av <- do.call(rbind, lapply(levels(cols), function(g) Rfast::colsums(cm[cols==g,]>0))) / (group.size+1);
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av;
  # fraction of cells expressing
  f.expressing <- t(gnzz / group.size);
  max.group <- max.col(log2.fold.change)

  if(verbose) {
    print("Summarizing results ... ")
  }

  ## sumarize
  ds <- lapply(1:ncol(x),function(i) {
    r <- data.frame(Z=x[,i],M=log2.fold.change[,i],highest=max.group==i,fe=f.expressing[,i])
    rownames(r) <- rownames(x)
    r
  })
  names(ds)<-colnames(x)

  return(ds)
}

## diffGenes is output of getDifferentialGenes
getMarkerGenes <- function(mat, com, diffGenes = NULL, upregulated.only=TRUE, z.threshold=3, auc.threshold=0.5, fe.threshold=0.5, M.threshold=0.1, verbose=TRUE, ncores=1) {

    if(is.null(diffGenes)) {
        ds <- getDifferentialGenes(mat, com)
    } else {
        ds <- diffGenes
    }

    if(upregulatedOnly) {
        dg <- unique(unlist(lapply(ds, function(x) rownames(x)[x$Z>z.threshold])))
    } else {
        dg <- unique(unlist(lapply(ds, function(x) rownames(x)[abs(x$Z)>z.threshold])))
    }
    dg <- intersect(dg, rownames(mat))
    mat <- mat[dg,]

    ## get info
    df.info <- df.list <- getClusterGeneInfo(mat, com, verbose, ncores)

    ## filter
    genes <- lapply(df.list, function(df) {
        if(upregulated.only) {
            df <- df[df$z.score > z.threshold,]
        } else {
            df <- df[abs(df$z.score) > z.threshold,]
        }
        df <- df[df$marker.auc > auc.threshold,]
        df <- df[df$fe > fe.threshold,]
        df <- df[df$M > M.threshold,]
        return(rownames(df))
    })
    names(genes) <- names(df.list)

    return(list(
        genes=genes,
        info=df.info
    ))
}

getStableClusters <- function(cm, cols, z.threshold=3, fc=2, t=0.6, verbose=TRUE, plot=TRUE) {
  if(verbose) {
    print(paste0('Feature selection for ', length(levels(cols)), ' groups with Z threshold ', z.threshold))
  }

  ## feature selection
  ds <- getDifferentialGenes(cm, cols, verbose=verbose)
  diff.genes <- lapply(ds, function(x) rownames(x)[abs(x$Z)>1.96 & x$highest])
  #lapply(diff.genes, length)
  ## visualize
  #lapply(diff.genes, function(dg) {
  #  dg <- intersect(dg, colnames(cm))
  #  m <- Rfast::rowsums(cm[,dg]>0)
  #  names(m) <- rownames(cm)
  #  plotEmbedding(myMudanObject1$emb[['PCA']], colors=m)
  #})

  diff.genes <- unique(unlist(diff.genes))
  diff.genes <- intersect(diff.genes, colnames(cm))
  length(diff.genes)

  ## assess accuracy
  set.seed(0)
  cells <- rownames(cm); folds <- split(cells, ceiling(seq_along(cells)/(length(cells)/fc)))
  #lapply(folds, length)

  ## subsample
  #set.seed(0)
  #train <- sample(rownames(cm), nrow(cm)*0.5)
  #test <- setdiff(rownames(cm), train)

  preds <- lapply(folds, function(train) {
    test <- setdiff(rownames(cm), train)
    df <- data.frame(celltype=cols, as.matrix(cm[, diff.genes]))
    ## model
    model <- MASS::lda(celltype ~ ., data=df[train,])
    ## predict
    model.output <- predict(model, df[test,])
    pred <- model.output$class
    names(pred) <- test
    return(pred)
  })
  preds <- unlist(preds)
  names(preds) <- sub('*.[.]','',names(preds))
  #head(preds)

  # ## assess accuracy
  # acc <- sapply(levels(cols), function(g) {
  #   #print(g)
  #   vi <- names(preds)[cols[names(preds)]==g]
  #   #print(table(preds[vi]))
  #   acc <- table(preds[vi]==cols[vi])
  #   acc <- acc[c('TRUE', 'FALSE')]
  #   names(acc) <- c('TRUE', 'FALSE')
  #   acc[is.na(acc)] <- 0
  #
  #   if(acc['TRUE']/sum(acc) < 0.5) {
  #     foo = table(preds[vi])
  #     foo = foo[-which(names(foo)==g)]
  #     print(paste0('May want to merge ', g , ' with ', names(sort(foo, decreasing=TRUE)[1])))
  #   }
  #
  #   return(acc['TRUE']/sum(acc))
  # })
  # names(acc) <- levels(cols)
  # print(acc)

  ## automerge
  newlevels <- sapply(levels(cols), function(g) {
    #print(g)
    vi <- names(preds)[cols[names(preds)]==g]
    #print(table(preds[vi]))
    acc <- table(preds[vi]==cols[vi])
    acc <- acc[c('TRUE', 'FALSE')]
    names(acc) <- c('TRUE', 'FALSE')
    acc[is.na(acc)] <- 0

    if(acc['TRUE']/sum(acc) < t) {
      foo = table(preds[vi])
      foo = foo[-which(names(foo)==g)]
      ng = names(sort(foo, decreasing=TRUE)[1])
      print(paste0('Merging ', g , ' with ', ng))
      return(ng)
    } else {
      return(g)
    }
  })

  cols2 <- factor(cols, levels=names(newlevels), labels=newlevels)
  cols2 <- factor(cols2)

  if(plot) {
    par(mfrow=c(1,2))
    plotEmbedding(myMudanObject1$emb[['PCA']], groups=cols, main='Original', show.legend=TRUE, mark.clusters=TRUE)
    plotEmbedding(myMudanObject1$emb[['PCA']], groups=cols2, main='Merged', show.legend=TRUE, mark.clusters=TRUE)
  }

  return(cols2)
}

##' Run tSNE on LDs from model
##'
##' @export
##'
tsneLda <- function(mat, model, perplexity=30, verbose=TRUE, plot=TRUE, do.par=TRUE, ncores=10, details=FALSE, ...) {
    if(verbose) {
        print('Running LDA models ...')
    }

    ## compute LDs
    if(class(model)=='lda') {
        ##reduction <- t(mat[rownames(model$scaling),]) %*% model$scaling
        reduction <- predict(model, data.frame(t(as.matrix(mat))))$x
    }
    if(class(model)=='list') {
        reduction <- do.call(cbind, lapply(model, function(m) {
                           reduction <- predict(m, data.frame(t(as.matrix(mat))))$x
                       }))
    }

    if(verbose) {
        print(paste0("Running Rtsne with perplexity ", perplexity))
    }

    ## tSNE
    emb <- Rtsne::Rtsne(reduction, is_distance=FALSE, perplexity=perplexity, verbose=verbose, num_threads=ncores)$Y
    rownames(emb) <- colnames(mat)

    if(plot) {
        if(do.par) {
            par(mar = c(0.5,0.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
        }
        plotEmbedding(emb, ...)
    }

    if(details) {
        return(list(
            emb=emb,
            reduction=reduction
        ))
    } else {
        return(emb)
    }
}


##' Plot 2D embedding
##'
##' @export
##'
plotEmbedding <- function(emb, groups=NULL, colors=NULL, do.par=TRUE, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, verbose=TRUE, unclassified.cell.color='gray70', group.level.colors=NULL, ...) {

    if(!is.null(colors)) {
        ## use clusters information
        if(!all(rownames(emb) %in% names(colors))) { warning("provided cluster vector doesn't list colors for all of the cells; unmatched cells will be shown in gray. ")}
        if(all(areColors(colors))) {
            if(verbose) cat("using supplied colors as is\n")
            cols <- colors[match(rownames(emb),names(colors))]; cols[is.na(cols)] <- unclassified.cell.color;
            names(cols) <- rownames(emb)
        } else {
            if(is.numeric(colors)) { # treat as a gradient
                if(verbose) cat("treating colors as a gradient")
                if(is.null(gradientPalette)) { # set up default gradients
                    if(all(sign(colors)>=0)) {
                        gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
                    } else {
                        gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
                    }
                }
                if(is.null(zlim)) { # set up value limits
                    if(all(sign(colors)>=0)) {
                        zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile)))
                        if(diff(zlim)==0) {
                            zlim <- as.numeric(range(colors))
                        }
                    } else {
                        zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile))
                        if(diff(zlim)==0) {
                            zlim <- c(-1,1)*as.numeric(max(abs(colors)))
                        }
                    }
                }
                                        # restrict the values
                colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];

                if(verbose) cat(' with zlim:',zlim,'\n')
                colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
                cols <- gradientPalette[colors[match(rownames(emb),names(colors))]*(length(gradientPalette)-1)+1]
                names(cols) <- rownames(emb)
            } else {
                stop("colors argument must be a cell-named vector of either character colors or numeric values to be mapped to a gradient")
            }
        }
    } else {
        if(!is.null(groups)) {
            if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- droplevels(groups); }
            groups <- as.factor(groups)[rownames(emb)]
            if(verbose) cat("using provided groups as a factor\n")
            factor.mapping=TRUE;
            ## set up a rainbow color on the factor
            factor.colors <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size,unclassified.cell.color=unclassified.cell.color,level.colors=group.level.colors,return.details=T)
            cols <- factor.colors$colors;
            names(cols) <- rownames(emb)
        } else {
          cols <- rep(unclassified.cell.color, nrow(emb))
          names(cols) <- rownames(emb)
        }
    }

    if(do.par) {
        par(mar = c(0.5,0.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    }
    plot(emb,col=adjustcolor(cols,alpha=alpha),cex=cex,pch=19,axes=F, ...); box();
    if(mark.clusters) {
        if(!is.null(groups)) {
            cent.pos <- do.call(rbind,tapply(1:nrow(emb),groups,function(ii) apply(emb[ii,,drop=F],2,median)))
            cent.pos <- na.omit(cent.pos);
            text(cent.pos[,1],cent.pos[,2],labels=rownames(cent.pos),cex=mark.cluster.cex)
        }
    }
    if(show.legend) {
        if(factor.mapping) {
            legend(x=legend.x,pch=rep(19,length(levels(groups))),bty='n',col=factor.colors$palette,legend=names(factor.colors$palette))
        }
    }
}
## a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.details=F,unclassified.cell.color='gray50',level.colors=NULL) {
    x <- as.factor(x);
    if(min.group.size>1) {
        x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
        x <- droplevels(x)
    }
    if(is.null(level.colors)) {
        col <- rainbow(length(levels(x)),s=s,v=v);
    } else {
        col <- level.colors[1:length(levels(x))];
    }
    names(col) <- levels(x);

    if(shuffle) col <- sample(col);

    y <- col[as.integer(x)]; names(y) <- names(x);
    y[is.na(y)] <- unclassified.cell.color;
    if(return.details) {
        return(list(colors=y,palette=col))
    } else {
        return(y);
    }
}
## quick utility to check if given character vector is colors
## thanks, Josh O'Brien: http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
    is.character(x) & sapply(x, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})
}







