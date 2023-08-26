#' Sample Methylation Data in Betas Format
#'
#' @description A subset of methylation data from The Cancer Genome Atlas (TCGA)
#' breast cancer cohort collected using 450k methylation arrays.
#'
#' @return A data frame with 293 rows and 300 columns, where each column
#' corresponds to a TCGA-BRCA sample and each row corresponds to an Illumina
#' methylation array probe.
#'
#' @details The dataset contains 293 rows representing individual probes and
#' 300 columns representing samples from the TCGA-BRCA cohort. The values in
#' the data frame represent methylation levels.
#'
#' @source \url{https://portal.gdc.cancer.gov/projects/TCGA-BRCA}
#' @usage data("betas")
#'
#' @examples
#' data(betas)
#' head(betas)
"betas"


#' Title: Annotations for Gene Promoter Regions
#'
#' @description Annotations containing the
#' start and end boundaries of gene promoters for protein-coding genes.
#'
#' @format ## `gene_annots`
#' A data frame with 19374 rows and 16 columns
#'
#' @return ## `gene_annots`
#' A data frame with 19374 rows and 16 columns
#' \describe{
#'   \item{seqnames}{Chromosome names (e.g., chr1, chr2, ...)}
#'   \item{start}{Start position of the promoter region}
#'   \item{end}{End position of the promoter region}
#'   \item{width}{Width of the promoter region}
#'   \item{strand}{Strand (+/-)}
#'   \item{tx_id}{Transcript ID (ENST)}
#'   \item{type}{Type of the region (e.g., hg38)}
#'   \item{gencode_gene_id}{Gencode Gene ID (ENSG)}
#'   \item{gencode_gene_type}{Gene type (e.g., protein coding)}
#'   \item{gencode_gene_name}{Gene name (e.g., SAMD11)}
#'   \item{transcript_type}{Type of the transcript (e.g., protein coding)}
#'   \item{transcript_name}{Transcript name}
#'   \item{transcript_support_level}{Transcript support level}
#'   \item{tag}{Tag for the entry}
#'   \item{is_canonical}{Logical indicating if the transcript is canonical}
#'   \item{gencode_region}{Gencode region information}
#' }
#'
#' @source None
#'
#' @usage data(gene_annots)
#'
#' @examples
#' data(gene_annots)
#' head(gene_annots)
"gene_annots"
