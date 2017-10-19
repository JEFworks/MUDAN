library(Rfast)
library(rbenchmark)

## generate random matrix
G=5; N=100; M=1000; initmean=0; initvar=10; upreg=10; upregvar=10; ng=100; seed=0
mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
mat[mat<0] <- 0
rownames(mat) <- paste0('gene', 1:M)
colnames(mat) <- paste0('cell', 1:(N*G))
##group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
group <- factor(sprintf('group%s', rep(1:G, each = N)))
names(group) <- colnames(mat)

## normalize
rbenchmark::benchmark(
                rowSums(mat),
                Rfast::rowsums(mat)
            )

mat <- t(t(mat)/Rfast::colsums(mat))
pseudocount <- 1
mat <- log10(mat*1e6+pseudocount)

## get variable genes
rbenchmark::benchmark(
            apply(mat, 1, var),
            Rfast::rowVars(mat)
            )

vi <- Rfast::rowVars(mat)
names(vi) <- rownames(mat)
n.odgenes <- 1000
odgenes <- names(sort(vi, decreasing=TRUE)[seq_len(n.odgenes)])

fast.pca <- function(m, nPcs=2, tol=1e-10, scale=FALSE, center=FALSE, transpose=FALSE) {
    require(irlba)
    ## note transpose is meant to speed up calculations when neither scaling nor centering is required
    if(transpose) {
        if(center) { m <- m-Matrix::rowMeans(m)}; if(scale) { m <- m/sqrt(Matrix::rowSums(m*m)); }
        a <- irlba(tcrossprod(m)/(ncol(m)-1), nu=0, nv=nPcs,tol=tol);
        a$l <- t(t(a$v) %*% m)
    } else {
        if(scale||center) { m <- scale(m,scale=scale,center=center) }
        a <- irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
        a$l <- m %*% a$v
    }
    a
}

nPcs <- 2
pcs <- fast.pca(t(mat[odgenes,]), nPcs=nPcs)
m <- t(pcs$l)
colnames(m) <- colnames(mat)
rownames(m) <- paste0('PC', seq_len(nPcs))

