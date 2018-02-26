#' Frozen PBMCs (Donor A)
#'
#' @format A sparse matrix with 13939 genes (rows) and 2896 cells (columns).
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"pbmcA"

#' Frozen PBMCs (Donor B)
#'
#' @format A sparse matrix with 15325 genes (rows) and 7765 cells (columns).
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"pbmcB"

#' Frozen PBMCs (Donor C)
#'
#' @format A sparse matrix with 16144 genes (rows) and 9493 cells (columns).
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"pbmcC"

#' Sorted CD14+ Monocytes CD19+ B Cells, CD34+ Cells, CD4+ Helper T Cells,
#' CD4+/CD25+ Regulatory T Cells, CD4+/CD45RA+/CD25- Naive T cells,
#' CD4+/CD45RO+ Memory T Cells, CD56+ Natural Killer Cells, CD8+ Cytotoxic
#' T cells, CD8+/CD45RA+ Naive Cytotoxic T Cells
#'
#' @format A sparse matrix with 8863 genes (rows) and 2140 cells (columns).
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"referenceCounts"

#' Annotations for sorted CD14+ Monocytes CD19+ B Cells, CD34+ Cells, CD4+
#' Helper T Cells, CD4+/CD25+ Regulatory T Cells, CD4+/CD45RA+/CD25- Naive
#' T cells, CD4+/CD45RO+ Memory T Cells, CD56+ Natural Killer Cells, CD8+
#' Cytotoxic T cells, CD8+/CD45RA+ Naive Cytotoxic T Cells
#'
#' @examples
#' \dontrun{
#' table(referenceAnnot)
#' #referenceAnnot
#' #        bcells     cytotoxict        memoryt      monocytes naivecytotoxic
#' #           330            148            252            129             43
#' #        naivet             nk    regulatoryt        thelper
#' #            52            769            226            191
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"referenceAnnot"
