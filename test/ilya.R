## exploring Ilya's idea of iterative refining
library(MUDAN)

## Analyze just PBMCA
data(pbmcA)
## restrict cells per sample to speed things up
A <- pbmcA[,1:1000]
## remove B cells from pbmcA identified by previous analysis
b.cells.pbmcA <- c("frozen_pbmc_donor_a_AAACATTGGCTAAC", "frozen_pbmc_donor_a_AAACCGTGTTACCT",
                   "frozen_pbmc_donor_a_AAACGCACTTCTAC", "frozen_pbmc_donor_a_AAAGACGAATGACC",
                   "frozen_pbmc_donor_a_AAAGGCCTCGAACT", "frozen_pbmc_donor_a_AAAGGCCTTGACCA",
                   "frozen_pbmc_donor_a_AAATCCCTAACCAC", "frozen_pbmc_donor_a_AAATGTTGAGATGA",
                   "frozen_pbmc_donor_a_AAATTCGAATTTCC", "frozen_pbmc_donor_a_AACACGTGCCCTTG",
                   "frozen_pbmc_donor_a_AACACTCTAGCAAA", "frozen_pbmc_donor_a_AACAGAGACTCGCT",
                   "frozen_pbmc_donor_a_AACCACGAAACAGA", "frozen_pbmc_donor_a_AACCACGAATGCTG",
                   "frozen_pbmc_donor_a_AACCGATGTGTGGT", "frozen_pbmc_donor_a_AACCGCCTGTTGTG",
                   "frozen_pbmc_donor_a_AACGCCCTTTCGGA", "frozen_pbmc_donor_a_AACGTGTGGTAGCT",
                   "frozen_pbmc_donor_a_AACTGTCTACCAAC", "frozen_pbmc_donor_a_AAGCAAGAGATAGA",
                   "frozen_pbmc_donor_a_AAGTAGGAAACCGT", "frozen_pbmc_donor_a_AATCCTTGCGAATC",
                   "frozen_pbmc_donor_a_AATCTCACCTCGAA", "frozen_pbmc_donor_a_AATGATACAGCACT",
                   "frozen_pbmc_donor_a_ACACCCTGAAAACG", "frozen_pbmc_donor_a_ACAGGTACTAGCGT",
                   "frozen_pbmc_donor_a_ACAGTCGAATGACC", "frozen_pbmc_donor_a_ACATTCTGGTACCA",
                   "frozen_pbmc_donor_a_ACCATTACTGTTTC", "frozen_pbmc_donor_a_ACCCAAGAGCCTTC",
                   "frozen_pbmc_donor_a_ACCTGGCTATCACG", "frozen_pbmc_donor_a_ACGAACACGAGCTT",
                   "frozen_pbmc_donor_a_ACGATTCTCCTACC", "frozen_pbmc_donor_a_ACGATTCTGTTCTT",
                   "frozen_pbmc_donor_a_ACGCACCTGTCATG", "frozen_pbmc_donor_a_ACGCACCTTTGGCA",
                   "frozen_pbmc_donor_a_ACGCTCACTGACAC", "frozen_pbmc_donor_a_ACGGATTGACGTAC",
                   "frozen_pbmc_donor_a_ACTATCACGGTGGA", "frozen_pbmc_donor_a_AGACACACCCCTTG",
                   "frozen_pbmc_donor_a_AGACTTCTCTGGAT", "frozen_pbmc_donor_a_AGATTAACCCTGAA",
                   "frozen_pbmc_donor_a_AGCATCGACCTTTA", "frozen_pbmc_donor_a_AGGACACTTAACCG",
                   "frozen_pbmc_donor_a_AGGGACGAAAACAG", "frozen_pbmc_donor_a_AGGTTGTGGAGGAC",
                   "frozen_pbmc_donor_a_AGTCACGAAAGATG", "frozen_pbmc_donor_a_AGTGTGACTTGACG",
                   "frozen_pbmc_donor_a_AGTTTAGATCGCCT", "frozen_pbmc_donor_a_AGTTTCACTTGCGA",
                   "frozen_pbmc_donor_a_ATAGGCTGAGTTCG", "frozen_pbmc_donor_a_ATAGGCTGGCGAGA",
                   "frozen_pbmc_donor_a_ATCAACCTGGTGGA", "frozen_pbmc_donor_a_ATCAGGTGGAGATA",
                   "frozen_pbmc_donor_a_ATCCAGGAGCAAGG", "frozen_pbmc_donor_a_ATCCCGTGTGTAGC",
                   "frozen_pbmc_donor_a_ATCCTAACCAGAAA", "frozen_pbmc_donor_a_ATCGCCACGAATCC",
                   "frozen_pbmc_donor_a_ATCTACTGCAGAAA", "frozen_pbmc_donor_a_ATCTCAACGGTTTG",
                   "frozen_pbmc_donor_a_ATGAGCACACCACA", "frozen_pbmc_donor_a_ATGGACACTGGAGG",
                   "frozen_pbmc_donor_a_ATGTACCTGTGTTG", "frozen_pbmc_donor_a_ATGTCACTGTCGAT",
                   "frozen_pbmc_donor_a_ATGTTGCTCGGGAA", "frozen_pbmc_donor_a_ATTAGTGAAGCGGA",
                   "frozen_pbmc_donor_a_ATTCAGCTGTCCTC", "frozen_pbmc_donor_a_ATTCTTCTGAAAGT",
                   "frozen_pbmc_donor_a_ATTTCCGATCAGTG", "frozen_pbmc_donor_a_CAAGGACTACCCAA",
                   "frozen_pbmc_donor_a_CAAGTTCTGTTTCT", "frozen_pbmc_donor_a_CACAGATGCCATAG",
                   "frozen_pbmc_donor_a_CACCACTGTGGTGT", "frozen_pbmc_donor_a_CACCGGGACTTGAG",
                   "frozen_pbmc_donor_a_CACCGTTGTCGTTT", "frozen_pbmc_donor_a_CACTCCGAGTCACA",
                   "frozen_pbmc_donor_a_CACTTATGTTGTGG", "frozen_pbmc_donor_a_CAGACTGAATAAGG",
                   "frozen_pbmc_donor_a_CAGACTGACGACTA")
A <- A[, setdiff(colnames(A), b.cells.pbmcA)]
dim(A)

## standard analysis
cd <- cleanCounts(A, min.reads = 10, min.detected = 10)
mat <- normalizeCounts(cd) ## CPM normalization
matnorm.info <- normalizeVariance(mat, details=TRUE) ## variance normalize
matnorm <- log10(matnorm.info$mat+1) ## log transform
pcs <- getPcs(matnorm[matnorm.info$ods,], nGenes=length(matnorm.info$ods), nPcs=30) ## 30 PCs on overdispersed genes
com <- getComMembership(pcs, k=10, method=igraph::cluster_infomap) ## graph-based community detection; over cluster with small k
emb <- Rtsne::Rtsne(pcs, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y ## get tSNE embedding on PCs
rownames(emb) <- rownames(pcs)
plotEmbedding(emb, com, mark.clusters=TRUE) ## plot
## get stable clusters
comA <- getStableClusters(cd, com, matnorm, min.group.size=10, min.diff.genes=10, z.threshold=1.96)
plotEmbedding(emb, groups=comA$com, mark.clusters=TRUE) ## plot
## train LDA model
genes <- intersect(unique(unlist(comA$pv.sig)), rownames(matnorm)) ## train on detected significantly differentially expressed genes among stable clusters
mn <- matnorm[genes,]
model <- modelLda(mat=mn, com=comA$com, retest=FALSE)
## remember normalization parameters
gsf <- matnorm.info$df$gsf
names(gsf) <- rownames(matnorm.info$df)
ods <- rownames(matnorm.info$mat)[matnorm.info$ods]
## project to LD space
preds <- predict(model, data.frame(t(log10(as.matrix(mat[names(gsf),]*gsf+1)))))
lds <- preds$x
class <- getConfidentPreds(preds$posterior)
emb.lds <- Rtsne::Rtsne(lds, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y ## get tSNE embedding on PCs
rownames(emb.lds) <- rownames(lds)
## final plot
par(mfrow=c(1,2))
plotEmbedding(emb, groups=class, mark.clusters=TRUE) ## plot
plotEmbedding(emb.lds, groups=class, mark.clusters=TRUE) ## plot
## save for later
A.info <- list(
  emb=emb,
  emb.lds=emb.lds,
  com=class,
  pv.sig=comA$pv.sig,
  gsf=gsf,
  ods=ods,
  mat=mat,
  model=model
)

## Bring in new sample; note this sample does have B cells
data(pbmcC)
C <- pbmcC[,1:2000]
dim(C)
## standard analysis
cd <- cleanCounts(C, min.reads = 10, min.detected = 10)
mat <- normalizeCounts(cd) ## CPM normalization
matnorm.info <- normalizeVariance(mat, details=TRUE) ## variance normalize
matnorm <- log10(matnorm.info$mat+1) ## log transform
pcs <- getPcs(matnorm[matnorm.info$ods,], nGenes=length(matnorm.info$ods), nPcs=30) ## 30 PCs on overdispersed genes
com <- getComMembership(pcs, k=30, method=igraph::cluster_infomap) ## graph-based community detection; over cluster with small k
emb <- Rtsne::Rtsne(pcs, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y ## get tSNE embedding on PCs
rownames(emb) <- colnames(matnorm)
plotEmbedding(emb, com, mark.clusters=TRUE) ## plot
## get stable clusters
comA <- getStableClusters(cd, com, matnorm, min.group.size=10, min.diff.genes=10, z.threshold=1.96)
graph.com <- comA$com
plotEmbedding(emb, groups=graph.com, mark.clusters=TRUE) ## plot

## predict A's annotations
predA <- predictLds(mat, A.info$model, A.info$gsf)
pred.com <- na.omit(getConfidentPreds(predA$posterior))
plotEmbedding(emb, groups=pred.com, mark.clusters=TRUE) ## plot

## combine A's predictions and new clusters
new.com <- factor(paste0(pred.com, '+', graph.com[names(pred.com)]))
names(new.com) <- names(pred.com)
table(new.com)
plotEmbedding(emb, groups=new.com, mark.clusters=TRUE, min.group.size=10) ## plot
com.comb <- getStableClusters(cd, new.com, matnorm, min.group.size=10, min.diff.genes=10, z.threshold=1.96)
plotEmbedding(emb, groups=as.factor(com.comb$com), mark.clusters=TRUE, min.group.size=10) ## plot

## train LDA model
genes <- intersect(unique(unlist(com.comb$pv.sig)), rownames(matnorm)) ## train on detected significantly differentially expressed genes among stable clusters
mn <- matnorm[genes,]
model <- modelLda(mat=mn, com=com.comb$com, retest=FALSE)
## remember normalization parameters
gsf <- matnorm.info$df$gsf
names(gsf) <- rownames(matnorm.info$df)
ods <- rownames(matnorm.info$mat)[matnorm.info$ods]
## project to LD space
preds <- predict(model, data.frame(t(log10(as.matrix(mat[names(gsf),]*gsf+1)))))
lds <- preds$x
class <- getConfidentPreds(preds$posterior)
emb.lds <- Rtsne::Rtsne(lds, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y ## get tSNE embedding on PCs
rownames(emb.lds) <- colnames(matnorm)
## final plot
par(mfrow=c(1,2))
plotEmbedding(emb, groups=class, mark.clusters=TRUE) ## plot
plotEmbedding(emb.lds, groups=class, mark.clusters=TRUE) ## plot
## save for later
C.info <- list(
  emb=emb,
  emb.lds=emb.lds,
  com=class,
  pv.sig=comA$pv.sig,
  gsf=gsf,
  ods=ods,
  mat=mat,
  model=model
)

## Use what we've learned from C to reanalyze A
predC <- predictLds(A.info$mat, C.info$model, C.info$gsf)
pred.comb <- getConfidentPreds(predC$posterior)
plotEmbedding(A.info$emb, groups=pred.comb, mark.clusters=TRUE) ## plot
## should not see any clusters predicted as B cells still

## get common embedding
cds <- list(A, C)
genes.int <- Reduce(intersect, lapply(cds, rownames))
cds.filtered <- lapply(cds, function(x) as.matrix(x[genes.int,]))
cds.all <- cleanCounts(do.call(cbind, cds.filtered), min.detected=10)
batch <- factor(unlist(lapply(colnames(cds.all), function(x) strsplit(x, '_')[[1]][4]))); names(batch) <- colnames(cds.all) ## get sample annotations
table(batch)
mat.all <- normalizeCounts(cds.all)
preds.all <- predictLds(mat.all, C.info$model, C.info$gsf)
lds.all <- preds.all$x
com.final <- factor(preds.all$class); names(com.final) <- rownames(preds.all$x)
## batch correct within clusters
lds.bc <- clusterBasedBatchCorrect(lds.all, batch, com.final)
## get embedding
emb.lds.bc <- Rtsne::Rtsne(lds.bc, is_distance=FALSE, perplexity=30, verbose=TRUE)$Y
rownames(emb.lds.bc) <- rownames(lds.bc)

## compare combined to individual embedding
par(mfrow=c(2,2), mar=c(0.5,0.5,2,0.5))
plotEmbedding(emb.lds.bc, com.final, alpha=0.1, main="MUDAN with combined annotations", show.legend=TRUE, legend.x = "topleft")
plotEmbedding(emb.lds.bc, batch, alpha=0.1, main="MUDAN with batch annotations", show.legend=TRUE, legend.x = "topleft")
plotEmbedding(A.info$emb, com.final, alpha=0.1, main="PBMCA", show.legend=TRUE, legend.x = "topleft")
plotEmbedding(C.info$emb, com.final, alpha=0.1, main="PBMCC", show.legend=TRUE, legend.x = "topleft")
