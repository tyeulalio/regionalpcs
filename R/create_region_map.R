#' Create a Region Map Between CpGs and Gene Regions
#'
#' This function generates a map that assigns CpG sites to gene regions,
#' establishing a linkage based on their genomic coordinates and providing
#' a foundation for subsequent region-specific analyses.
#'
#' @param cpg_gr A `GRanges` object containing the genomic positions of CpG
#' sites.
#' @param genes_gr A `GRanges` object containing the genomic positions of gene
#' regions (e.g., promoters) of interest.
#' @param verbose Boolean; print output statements
#'
#' @return A `data.frame` with mappings between gene IDs and CpG IDs,
#' facilitating associating CpG sites with their corresponding gene regions for
#'  downstream analyses.
#'
#' @export
#'
#' @examples
#' library(GenomicRanges)
#'
#' # Creating dummy GRanges objects for CpG sites and gene regions
#' cpg_gr <- GRanges(seqnames=c("chr1", "chr1", "chr2"),
#'                     ranges=IRanges(start=c(100, 200, 150),
#'                     end=c(100, 200, 150)))
#' genes_gr <- GRanges(seqnames=c("chr1", "chr2", "chr2"),
#'                     ranges=IRanges(start=c(50, 100, 130),
#'                     end=c(150, 180, 160)))
#
#' # Creating a region map using the function
#' region_map <- create_region_map(cpg_gr, genes_gr)
create_region_map <- function(cpg_gr, genes_gr, verbose=FALSE){
    # Check if GRanges objects have names and assign if absent
    if (is.null(names(cpg_gr))) {
        warning("'cpg_gr' does not have names, assigning default names.")
        names(cpg_gr) <- paste0("cpg", seq_along(cpg_gr))
    }
    if (is.null(names(genes_gr))) {
        warning("'genes_gr' does not have names, assigning default names.")
        names(genes_gr) <- paste0("gene", seq_along(genes_gr))
    }

    # Identify overlaps between CpG positions and gene regions
    # Utilizing the 'findOverlaps' function, the relationship between
    # CpGs and their corresponding gene regions is discerned,
    # based on their genomic coordinates.
    overlaps <- GenomicRanges::findOverlaps(query=genes_gr, subject=cpg_gr) |>
        as.data.frame()

    # Construct a region map
    # The dataframe 'region_map' is crafted to establish a correspondence
    # between gene IDs and CpG IDs, anchored in their genomic co-localization,
    # serving as a key for correlating methylation data with gene regions
    # in downstream analyses.
    region_map <- data.frame(gene_id=names(genes_gr)[overlaps$queryHits],
                                cpg_id=names(cpg_gr)[overlaps$subjectHits])

    # Optionally print the first few lines of the region_map
    if (verbose == TRUE) {
        message("Preview of the region map:")
        message(utils::head(region_map))
    }

    return(region_map)
}
