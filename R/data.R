#' Sample methylation data in betas format
#'
#' A subset of methylation data from The Cancer Genome Atlas (TCGA)
#' breast cancer cohort collected using 450k methylation arrays
#' Report ...
#'
#' @format ## `betas`
#' A data frame with 293 rows and 300 columns:
#' \describe{
#'   \item{rows}{Illumina methylation array probes}
#'   \item{columns}{TCGA-BRCA samples}
#'   ...
#' }
#' @source <https://portal.gdc.cancer.gov/projects/TCGA-BRCA>
"betas"


#' Annotations for gene promoter regions
#'
#' Annotations containing the start and end boundaries of
#' gene promoters for protein coding genes.
#' Report ...
#'
#' @format ## `gene_annots`
#' A data frame with XXX rows and XXX columns:
#' \describe{
#'   \item{country}{Country name}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{Year}
#'   ...
#' }
#' @source <>
"gene_annots"
