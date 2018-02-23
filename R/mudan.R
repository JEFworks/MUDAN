#' @title Filter a counts matrix
#' @description Filter a counts matrix based on gene (row) and cell (column)
#'      requirements.
#' @param counts A sparse read count matrix. The rows correspond to genes,
#'      columns correspond to individual cells
#' @param min.lib.size Minimum number of genes detected in a cell. Cells with
#'      fewer genes will be removed (default: 300)
#' @param max.lib.size Maximum number of genes detected in a cell. Cells with
#'      more genes will be removed (default: 8000)
#' @param min.reads Minimum number of reads per gene. Genes with fewer reads
#'      will be removed (default: 1)
#' @param min.detected Minimum number of cells a gene must be seen in. Genes
#'      not seen in a sufficient number of cells will be removed (default: 1)
#' @param verbose Verbosity (default: TRUE)
#' @return a filtered read count matrix
#' @examples
#' data(pbmcA)
#' dim(pbmcA)
#' mat <- cleanCounts(pbmcA)
#' dim(mat)
#' @export
#' @importFrom Matrix Matrix colSums rowSums
cleanCounts <- function(
  counts,
  min.lib.size = 300,
  max.lib.size = 8000,
  min.reads    = 1,
  min.detected = 1,
  verbose      = TRUE
) {

  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix")) {
    message('Converting to sparse matrix ...')
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  if (verbose) {
    message('Filtering matrix with ', ncol(counts), ' cells and ', nrow(counts), ' genes ...')
  }

  ## filter out low-gene cells or potential doublets
  ix_col <- Matrix::colSums(counts)
  ix_col <- ix_col > min.lib.size & ix_col < max.lib.size
  counts <- counts[, ix_col]

  ## remove genes that don't have many reads
  counts <- counts[Matrix::rowSums(counts) > min.reads, ]

  ## remove genes that are not seen in a sufficient number of cells
  counts <- counts[Matrix::rowSums(counts > 0) > min.detected, ]

  if (verbose) {
    message('Resulting matrix has ', ncol(counts), ' cells and ', nrow(counts), ' genes')
  }

  return(counts)
}


#' Normalizes counts to CPM
#'
#' Normalizes raw counts to log10 counts per million with pseudocount
#'
#' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param depthScale Depth scaling. Using a million for CPM (default: e61)
#' @param verbose Verbosity (default: TRUE)
#'
#' @return a normalized matrix
#'
#' @examples {
#' data(pbmcA)
#' cd <- pbmcA[, 1:500]
#' mat <- normalizeCounts(cd)
#' }
#'
#' @export
#'
normalizeCounts <- function(counts, depthScale=1e6, verbose=TRUE) {

  if(verbose) {
    print(paste0('Normalizing matrix with ', ncol(counts), ' cells and ', nrow(counts), ' genes'))
  }

  counts <- t(t(counts)/colSums(counts))
  counts <- counts*depthScale

  return(counts)
}


#' Normalize gene expression variance relative to transcriptome-wide expectations
#' (Modified from SCDE/PAGODA2 code)
#'
#' Normalizes gene expression magnitudes to with respect to its ratio to the
#' transcriptome-wide expectation as determined by local regression on expression magnitude
#'
#' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param gam.k Generalized additive model parameter; the dimension of the basis used to represent the smooth term (default: 5)
#' @param alpha Significance threshold (default: 0.05)
#' @param plot Whether to plot the results (default: FALSE)
#' @param use.unadjusted.pvals If true, will apply BH correction (default: FALSE)
#' @param do.par Whether to adjust par for plotting if plotting (default: TRUE)
#' @param max.adjusted.variance Ceiling on maximum variance after normalization to prevent infinites (default: 1e3)
#' @param min.adjusted.variance Floor on minimum variance after normalization (default: 1e-3)
#' @param verbose Verbosity (default: TRUE)
#' @param details If true, will return data frame of normalization parameters. Else will return normalized matrix.(default: FALSE)
#'
#' @return If details is true, will return data frame of normalization parameters. Else will return normalized matrix.
#'
#' @examples {
#' data(pbmcA)
#' cd <- pbmcA[, 1:500]
#' mat <- cleanCounts(cd)
#' mat <- normalizeVariance(mat)
#' }
#'
#' @importFrom mgcv s
#'
#' @export
#'
normalizeVariance <- function(counts, gam.k=5, alpha=0.05, plot=FALSE, use.unadjusted.pvals=FALSE, do.par=TRUE, max.adjusted.variance=1e3, min.adjusted.variance=1e-3, verbose=TRUE, details=FALSE) {
  mat <- t(counts) ## make rows as cells, cols as genes

  if(verbose) {
    print("Calculating variance fit ...")
  }
  dfm <- log(colMeans(mat))
  dfv <- log(apply(mat, 2, var))
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

  ## variance normalize
  norm.mat <- counts*df$gsf
  if(!details) {
    return(norm.mat)
  } else {
    ## return normalization factor
    return(list(mat=norm.mat, ods=ods, df=df))
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


#' Dimensionality reduction by PCA
#'
#' Dimensionality reduction using PCA by computing principal components using highly variable genes
#'
#' @param mat Variance normalized gene expression matrix.
#' @param nGenes Number of most variable genes. (default: 1000)
#' @param nPcs Number of principal components. (default: 100)
#' @param verbose Verbosity (default: TRUE)
#' @param ... Additional parameters to pass to irlba
#'
#' @return Matrix with columns as cells and rows as principal component eigenvectors.
#'
#' @examples {
#' data(pbmcA)
#' cd <- pbmcA[, 1:500]
#' mat <- cleanCounts(cd)
#' mat <- normalizeVariance(mat)
#' pcs <- getPcs(mat)
#' }
#'
#' @export
#'
getPcs <- function(mat, nGenes = min(nrow(mat), 1000), nPcs = 100, verbose=TRUE, ...) {
  if(class(mat)!='matrix') {
    mat <- as.matrix(mat)
  }

  if(verbose) {
    print(paste0('Identifying top ', nPcs, ' PCs on ', nGenes, ' most variable genes ...'))
  }

  ## get variable genes
  if(nGenes < nrow(mat)) {
    vgenes <- getVariableGenes(mat, nGenes)
    mat <- mat[vgenes,]
  }

  ## get PCs
  pcs <- fastPca(t(mat), nPcs, center=TRUE, ...)
  m <- t(pcs$l)
  colnames(m) <- colnames(mat)
  rownames(m) <- paste0('PC', seq_len(nPcs))

  return(t(m))
}
#' Derive principal components by SVD
#'
#' @param m Matrix X
#' @param nPcs Number of principal components
#' @param tol irlba tol
#' @param scale Scale input
#' @param center Center input
#' @param ... Additional parameters to pass to irlba
#'
#' @export
#'
fastPca <- function(m, nPcs=2, tol=1e-10, scale=FALSE, center=TRUE, ...) {
  if(scale||center) {
    m <- scale(m, scale=scale, center=center)
  }
  a <- irlba::irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs, tol=tol, ...)
  a$l <- m %*% a$v

  return(a)
}
#' Get variable genes
#'
#' @param mat Gene expression matrix
#' @param nGenes Number of genes
#'
#' @export
#'
getVariableGenes <- function(mat, nGenes) {
  vi <- apply(mat, 1, var)
  names(vi) <- rownames(mat)
  vgenes <- names(sort(vi, decreasing=TRUE)[seq_len(nGenes)])
  return(vgenes)
}


#' Get clusters by community detection on approximate nearest neighbors
#'
#' Group cells into clusters based on graph-based community detection on approximate nearest neighbors
#'
#' @param mat Matrix of cells as columns. Features as rows (such as PCs).
#' @param k K-nearest neighbor parameter.
#' @param method Community detection method from igraph. (default: igraph::cluster_walktrap)
#' @param verbose Verbosity (default: TRUE)
#' @param details Whether to return just community annotation or entended details including graph and graph modularity
#'
#' @return Vector of community annotations
#'
#' @examples {
#' data(pbmcA)
#' cd <- pbmcA[, 1:500]
#' mat <- cleanCounts(cd)
#' mat <- normalizeVariance(mat)
#' pcs <- getPcs(mat)
#' com <- getComMembership(pcs, k=30)
#' }
#'
#' @export
#'
getComMembership <- function(mat, k, method=igraph::cluster_walktrap, verbose=TRUE, details=FALSE) {
  if(verbose) {
    print("finding approximate nearest neighbors ...")
  }
  knn <- RANN::nn2(mat, k=k)[[1]]

  ## convert to adjacency matrix
  adj <- matrix(0, nrow(mat), nrow(mat))
  rownames(adj) <- colnames(adj) <- rownames(mat)
  invisible(lapply(seq_len(nrow(mat)), function(i) {
    adj[i,rownames(mat)[knn[i,]]] <<- 1
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


#' Get clusters by community detection on approximate nearest neighbors with subsampling
#'
#' Group cells into clusters based on graph-based community detection on approximate nearest neighbors for random subset of cells
#' For when getComMembership takes too long due to there being too many cells
#'
#' @param mat Matrix of cells as columns. Features as rows (such as PCs).
#' @param k K-nearest neighbor parameter.
#' @param nsubsample Number of cells in subset (default: ncol(mat)*0.5)
#' @param seed Random seed for reproducibility
#' @param vote Use neighbor voting system to annotate rest of cells not in subset. If false, will use machine-learning model. (default: FALSE)
#' @param method Community detection method from igraph. (default: igraph::cluster_walktrap)
#' @param verbose Verbosity (default: TRUE)
#'
#' @return Vector of community annotations
#'
#' @examples \dontrun{
#' data(pbmcA)
#' cd <- pbmcA
#' mat <- cleanCounts(cd)
#' mat <- normalizeVariance(mat)
#' pcs <- getPcs(mat)
#' com <- getApproxComMembership(pcs, k=30, getApproxComMembership=1000)
#' }
#'
#' @export
#'
getApproxComMembership <- function(mat, k, nsubsample=ncol(mat)*0.5, method=igraph::cluster_walktrap, seed=0, vote=FALSE, verbose=TRUE) {

  if(verbose) {
    print(paste0('Subsampling from ', ncol(mat), ' cells to ', nsubsample, ' ... '))
  }
  ## random subsampling
  ## TODO: density based downsampling
  set.seed(seed)
  subsample <- sample(rownames(mat), nsubsample)

  if(verbose) {
    print('Identifying cluster membership for subsample ... ')
  }
  pcs.sub <- mat[subsample, ]
  com.sub <- getComMembership(pcs.sub, k=k, method=method)

  if(verbose) {
    print('Imputing cluster membership for rest of cells ... ')
  }

  ## Use neighbor voting
  if(vote) {
    data <- mat[subsample, ]
    query <- mat[setdiff(rownames(mat), subsample), ]
    knn <- RANN::nn2(data, query, k=k)[[1]]
    rownames(knn) <- rownames(query)
    com.nonsub <- unlist(apply(knn, 1, function(x) {
      ## nearest neighbors in data
      nn <- rownames(data)[x]
      ## look at their cell type annotations
      nn.com <- com.sub[nn]
      ## get most frequent annotation
      return(names(sort(table(nn.com), decreasing=TRUE)[1]))
    }))
    com.all <- factor(c(com.sub, com.nonsub)[rownames(mat)])
  }
  else {
    ## Use model instead
    ## Inspired by DenSVM
    df.sub <- data.frame(celltype=com.sub, pcs.sub)
    model <- MASS::lda(celltype ~ ., data=df.sub)
    df.all <- data.frame(mat)
    model.output <- predict(model, df.all)
    com.all <- model.output$class
    names(com.all) <- rownames(df.all)
    if(verbose) {
      print("Model accuracy for subsample ...")
      print(table(com.all[names(com.sub)]==com.sub))
    }
  }

  return(com.all)
}


#' Linear discriminant analysis model
#'
#' Identifies components that maximally discriminate among groups using a linear discriminant analysis model
#'
#' @param mat Expression matrix with cells as columns, transferable features such as genes as rows.
#' @param com Community annotations
#' @param verbose Verbosity (default: TRUE)
#' @param nfeatures Number of features (genes) in LDA model (default: all)
#' @param random Wehther those features are random of chosen based on variance (most variable will be chosen by default)
#' @param retest Whether to retest model for accuracy
#'
#' @return LDA model
#'
#' @examples {
#' data(pbmcA)
#' cd <- pbmcA[, 1:500]
#' mat <- cleanCounts(cd)
#' mat <- normalizeVariance(mat)
#' pcs <- getPcs(mat)
#' com <- getComMembership(pcs, k=30)
#' model <- modelLda(mat, com)
#' }
#'
#' @export
#'
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


#' Differential expression analysis (adapted from PAGODA2)
#'
#' @param cd A read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param cols Column/cell group annotations. Will perform one vs. all differential expression analysis.
#' @param verbose Verbosity
#'
#' @export
#'
getDifferentialGenes <- function(cd, cols, verbose=TRUE) {
  cm <- t(cd)

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
    bar <- rank(foo)
    #bar[is.na(foo)] <- 0
    bar[foo==0] <- 0
    bar
  }); rownames(xr) <- rownames(cm)
  ##xr <- sparse_matrix_column_ranks(cm);

  # calculate rank sums per group
  grs <- do.call(rbind, lapply(levels(cols), function(g) colSums(xr[cols==g,])))
  rownames(grs) <- levels(cols); colnames(grs) <- colnames(xr)
  ##grs <- colSumByFac(xr,as.integer(cols))[-1,,drop=F]

  # calculate number of non-zero entries per group
  gnzz <- do.call(rbind, lapply(levels(cols), function(g) colSums(xr[cols==g,]>0)))
  rownames(gnzz) <- levels(cols); colnames(gnzz) <- colnames(xr)
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
    bar <- rank(foo)
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
  log.gene.av <- log2(colMeans(cm));
  group.gene.av <- do.call(rbind, lapply(levels(cols), function(g) colSums(cm[cols==g,]>0))) / (group.size+1);
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av;
  # fraction of cells expressing
  f.expressing <- t(gnzz / group.size);
  max.group <- max.col(log2.fold.change)

  if(verbose) {
    print("Summarizing results ... ")
  }

  ## summarize
  ds <- lapply(1:ncol(x),function(i) {
    r <- data.frame(Z=x[,i],M=log2.fold.change[,i],highest=max.group==i,fe=f.expressing[,i])
    rownames(r) <- rownames(x)
    r
  })
  names(ds)<-colnames(x)

  return(ds)
}
## annotate clusters; old; needs major speed improvement
getClusterGeneInfo <- function(mat, com, verbose=FALSE, ncores=10) {

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
    xm <- rowMeans(x)
    ym <- rowMeans(y)
    fc <- log2(xm/ym)
    ## percent expressing
    pe <- rowSums(x>0)/ncol(x)
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
## too slow; to do
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
## diffGenes is output of getDifferentialGenes
getMarkerGenes <- function(mat, com, diffGenes = NULL, upregulated.only=TRUE, z.threshold=3, auc.threshold=0.5, fe.threshold=0.5, M.threshold=0.1, verbose=TRUE, ncores=1) {

  if(is.null(diffGenes)) {
    ds <- getDifferentialGenes(mat, com)
  } else {
    ds <- diffGenes
  }

  if(upregulated.only) {
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


#' Iterative merging of clusters along tree until stable
#'
#' @param cd Counts matrix for differential expression analysis. Rows are genes. Columns are cells.
#' @param com Community/group annotations for cells
#' @param matnorm Normalized gene expression matrix for building group relationship tree. Rows are genes. Columns are cells.
#' @param z.threshold Z-score threshold for identifying significantly differentially expressed genes
#' @param hclust.method Hierarchical clustering method used to construct relationship tree
#' @param min.group.size Minimum group size for stable cluster
#' @param min.diff.genes Minimum number of significantly differentially expressed genes that must be identified or else groups will be merged
#' @param max.iter Maximum number of iterations. Will end earlier if convergence reached.
#' @param plot Whether to plot intermediate plots (hierarchical clustering dendrograms)
#' @param verbose Verbosity
#'
#' @export
#'
getStableClusters <- function(cd, com, matnorm, z.threshold=3, hclust.method='ward.D', min.group.size=10, min.diff.genes=nrow(cd)*0.005, max.iter=10, verbose=TRUE, plot=FALSE) {

  if(min.group.size>1) { com[com %in% levels(com)[unlist(tapply(com,com,length))<min.group.size]] <- NA; com <- droplevels(com); }
  com <- as.factor(com)[colnames(cd)]

  ## compute differentially expressed genes between two branches of dendrogram
  compare <- function(dend) {
    g1 <- labels(dend[[1]])
    g2 <- labels(dend[[2]])
    incom <- names(com)[com %in% g1]
    outcom <- names(com)[com %in% g2]
    com.sub <- factor(c(
      rep(paste0(g1, collapse="."), length(incom)),
      rep(paste0(g2, collapse="."), length(outcom))))
    names(com.sub) <- c(incom, outcom)
    dg <- getDifferentialGenes(cd[, names(com.sub)], com.sub, verbose=verbose)
    dg.sig <- lapply(dg, function(x) {
      ## get only significantly upregulated
      x <- x[x$Z>z.threshold,]
      x <- na.omit(x)
      return(rownames(x))
    })
    #dg.sig <- dg.sig[!grepl('RP|MT', dg.sig)] ## exclude ribosomal and mitochondrial?
    return(dg.sig)
  }
  ## recursively trace dendrogram
  pv.recur <- function(dend) {
    ## compares leaves of tree (groups)
    if(is.leaf(dend[[1]]) & is.leaf(dend[[2]])) {
      g1 <- paste0(labels(dend[[1]]), collapse=":")
      g2 <- paste0(labels(dend[[2]]), collapse=":")
      if(verbose) {
        print(paste0('Comparing ', g1, ' and ', g2))
      }
      ## if insufficient number of marker genes distinguishing two groups, merge
      dg.sig <- compare(dend)
      if(verbose) {
        print('Differential genes found: ')
        print(sapply(dg.sig, length))
      }
      ## both leaves need markers
      #if(any(sapply(dg.sig, length) < min.diff.genes)) {
      if(length(unlist(dg.sig)) < min.diff.genes) {
        if(verbose) {
          print(paste0('Merging ', g1, ' and ', g2))
        }
        com.fin[com.fin==g1] <<- paste0(g1,'.', g2)
        com.fin[com.fin==g2] <<- paste0(g1,'.', g2)
      }
      pv.sig.all[[g1]] <<- dg.sig
      pv.sig.all[[g2]] <<- dg.sig
    } else {
      ## if only one is leaf, compare to entire other branch
      if(is.leaf(dend[[1]])) {
        g1 <- paste0(labels(dend[[1]]), collapse=".")
        g2 <- paste0(labels(dend[[2]]), collapse=".")
        if(verbose) {
          print(paste0('Comparing ', g1, ' and ', g2))
        }
        dg.sig <- compare(dend)
        if(verbose) {
          print('Differential genes found: ')
          print(sapply(dg.sig, length))
        }
        ## if insufficient number of marker genes for leaf, set to NA
        #if(any(sapply(dg.sig, length) < min.diff.genes)) {
        #if(length(unlist(dg.sig)) < min.diff.genes) {
        if(length(dg.sig[[labels(dend[[1]])]]) < min.diff.genes) {
          if(verbose) {
            print(paste0('Cannot distinguish ', g1, ' and ', g2))
          }
          com.fin[com.fin==g1] <<- NA
        }
        pv.sig.all[[g1]] <<- dg.sig
      } else {
        pv.recur(dend[[1]])
      }
      if(is.leaf(dend[[2]])) {
        g1 <- paste0(labels(dend[[1]]), collapse=".")
        g2 <- paste0(labels(dend[[2]]), collapse=".")
        if(verbose) {
          print(paste0('Comparing ', g1, ' and ', g2))
        }
        dg.sig <- compare(dend)
        if(verbose) {
          print('Differential genes found: ')
          print(sapply(dg.sig, length))
        }
        #if(any(sapply(dg.sig, length) < min.diff.genes)) {
        #if(length(unlist(dg.sig)) < min.diff.genes) {
        if(length(dg.sig[[labels(dend[[2]])]]) < min.diff.genes) {
          if(verbose) {
            print(paste0('Cannot distinguish ', g1, ' and ', g2))
          }
          com.fin[com.fin==g2] <<- NA
        }
        pv.sig.all[[g2]] <<- dg.sig
      } else {
        pv.recur(dend[[2]])
      }
    }
  }

  i <- 1 ## counter
  repeat {
    com <- factor(com)
    ## average within groups
    mat.summary <- do.call(cbind, lapply(levels(com), function(ct) {
      cells <- which(com==ct)
      rowMeans(matnorm[, cells])
    }))
    colnames(mat.summary) <- levels(com)
    dim(mat.summary)
    ## cluster groups
    hc <- hclust(dist(t(mat.summary)), method=hclust.method)
    if(plot) { plot(hc) }
    dend <- as.dendrogram(hc)
    com.fin <- as.character(com)
    names(com.fin) <- names(com)
    pv.sig.all <- list()
    pv.recur(dend)
    ## test if converged
    if(i>max.iter) {
      break
    }
    if(sum(com==com.fin, na.rm=TRUE)>=min(sum(!is.na(com)), sum(!is.na(com.fin)))) {
      break
    }
    com <- com.fin
    i <- i+1
  }

  return(list(com=com.fin, pv.sig=pv.sig.all, hc=hc, mat.summary=mat.summary))

}


#' Predict LD embedding for new dataset given old model and gene scale factor
#'
#' @param mat Library-size normalized gene expression matrix
#' @param model LDA model
#' @param gsf Gene scale factor to be applied to mat (so mat must not be variance normalized)
#' @param verbose Verbosity
#'
#' @export
#'
predictLds <- function(mat, model, gsf, verbose=TRUE) {
  gsf.have <- intersect(names(gsf), rownames(mat))
  if(verbose) {
    print(paste0('Percentage of features retained: ', length(gsf.have)/length(gsf)))
  }

  ## fill in with zeros
  mat.test <- matrix(0, length(gsf), ncol(mat))
  rownames(mat.test) <- names(gsf)
  colnames(mat.test) <- colnames(mat)
  mat.test[gsf.have,] <- as.matrix(mat[gsf.have,])

  pred <- predict(model, data.frame(t(log10(mat.test*gsf+1))))
}


#' Use LDA model posteriors to retain only confident predictions
#'
#' @param posterior Posterior probabilities from LDA
#' @param t Posteriors below this threshold are set to NA
#'
#' @export
#'
getConfidentPreds <- function(posterior, t=0.95) {
  class <- apply(posterior, 1, function(x) {
    x <- sort(x, decreasing=TRUE)
    if(x[1]-sum(x[-1]) >= t) {
      return(names(x)[1])
    } else {
      return(NA)
    }
  })
  class <- factor(class)
  return(class)
}


#' Batch correct within identified groups using ComBat
#'
#' @param lds.all Matrix to be batch corrected
#' @param batch Batch factor annotations
#' @param com.final Group annotations
#' @param min.group.size Minimum number of cells in a group in order to batch correct
#' @export
#'
clusterBasedBatchCorrect <- function(lds.all, batch, com.final, min.group.size=10) {
  com.final <- factor(com.final)
  lds.bc <- do.call(rbind, lapply(levels(com.final), function(ct){
    cells <- na.omit(names(com.final)[com.final==ct])
    batch.cells <- factor(batch[cells])
    ## correct only if more than N cells per group
    if(sum(table(batch.cells)>min.group.size)>1) {
      t(sva::ComBat(t(lds.all[cells,]), batch.cells))
    } else {
      lds.all[cells,]
    }
  }))
  return(lds.bc)
}


#' Run tSNE on LDs from model
#'
#' @param mat Normalized matrix
#' @param model LDA model
#' @param perplexity Perplexity parameter for tSNE
#' @param verbose Verbosity
#' @param plot Whether to plot
#' @param do.par Whether to set plot margins
#' @param ncores Number of cores for paralele tSNE
#' @param details Whether to return details
#' @param ... Additional parameters to pass to plot
#'
#' @export
#'
tsneLda <- function(mat, model, perplexity=30, verbose=TRUE, plot=TRUE, do.par=TRUE, ncores=10, details=FALSE, ...) {
  if(verbose) {
    print('Running LDA models ...')
  }

  ## compute LDs
  matm <- t(as.matrix(mat))
  if(class(model)=='lda') {
    genes.need <- rownames(model$scaling)
    mat.temp <- matrix(0, nrow(matm), length(genes.need))
    rownames(mat.temp) <- rownames(matm)
    colnames(mat.temp) <- genes.need
    genes.have <- intersect(genes.need, colnames(matm))
    mat.temp[rownames(matm), genes.have] <- matm[, genes.have]
    predict(model, data.frame(mat.temp))$x
  }
  if(class(model)=='list') {
    reduction <- do.call(cbind, lapply(model, function(m) {
      genes.need <- rownames(m$scaling)
      mat.temp <- matrix(0, nrow(matm), length(genes.need))
      rownames(mat.temp) <- rownames(matm)
      colnames(mat.temp) <- genes.need
      genes.have <- intersect(genes.need, colnames(matm))
      mat.temp[rownames(matm), genes.have] <- matm[, genes.have]
      predict(m, data.frame(mat.temp))$x
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


#' Plot 2D embedding
#'
#' @param emb dataframe with x and y coordinates
#' @param groups factor annotations for rows on emb for visualizing cluster annotations
#' @param colors color or numeric values for rows on emb for visualizing gene expression
#' @param cex point size
#' @param alpha point opacity
#' @param gradientPalette palette for colors if numeric values provided
#' @param zlim range for colors
#' @param s saturation of rainbow for group colors
#' @param v value of rainbow for group colors
#' @param min.group.size minimum size of group in order for group to be colored
#' @param show.legend whether to show legend
#' @param mark.clusters whether to mark clusters with name of cluster
#' @param mark.cluster.cex cluster marker point size
#' @param shuffle.colors whether to shuffle group colors
#' @param legend.x legend position ie. 'topright', 'topleft', 'bottomleft', 'bottomright'
#' @param gradient.range.quantile quantile for mapping colors to gradient palette
#' @param verbose verbosity
#' @param unclassified.cell.color cells not included in groups will be labeled in this color
#' @param group.level.colors set group level colors. Default uses rainbow.
#' @param ... Additional parameters to pass to BASE::plot
#'
#' @export
#'
plotEmbedding <- function(emb, groups=NULL, colors=NULL, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, verbose=TRUE, unclassified.cell.color='gray70', group.level.colors=NULL, ...) {

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







