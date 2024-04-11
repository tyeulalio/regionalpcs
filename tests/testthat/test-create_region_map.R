# Creating dummy GRanges objects for CpG sites and gene regions
cpg_gr <- GenomicRanges::GRanges(seqnames=c("chr1", "chr1", "chr2"),
                                    ranges=IRanges::IRanges(start=c(100, 200, 150),
                                    end=c(101, 201, 151)))
genes_gr <- GenomicRanges::GRanges(seqnames=c("chr1", "chr2", "chr2"),
                                    ranges=IRanges::IRanges(start=c(50, 100, 130),
                                    end=c(150, 180, 160)))
names(cpg_gr) <- c('cpg1', 'cpg2', 'cpg3')
names(genes_gr) <- c('gene1', 'gene2', 'gene3')

test_that("create_region_map returns the expected output", {
    # Checking the basic function output type and dimensions
    region_map <- create_region_map(cpg_gr, genes_gr, verbose=TRUE)
    expect_s3_class(region_map, "data.frame")
    expect_named(region_map, c("gene_id", "cpg_id"))
    expect_equal(nrow(region_map), 3)

    # Verifying the content of the map based on genomic overlap logic
    expect_equal(region_map$gene_id, c("gene1", "gene2", "gene3"))
    expect_equal(region_map$cpg_id, c("cpg1", "cpg3", "cpg3"))
})

test_that("create_region_map manages no overlap scenario gracefully", {
    # Creating GRanges objects that have no overlapping regions
    cpg_gr_no_overlap <- GenomicRanges::GRanges(
        seqnames=c("chr1", "chr1", "chr2"),
        ranges=IRanges::IRanges(start=c(300, 400, 350),
        end=c(300, 400, 350)))
    genes_gr_no_overlap <- GenomicRanges::GRanges(
        seqnames=c("chr1", "chr2", "chr2"),
        ranges=IRanges::IRanges(start=c(50, 100, 130),
        end=c(150, 180, 160)))
    names(cpg_gr_no_overlap) <- c('cpg1', 'cpg2', 'cpg3')
    names(genes_gr_no_overlap) <- c('gene1', 'gene2', 'gene3')

    # Ensuring the function handles non-overlapping regions without
    # errors and gives empty df
    region_map_no_overlap <- create_region_map(cpg_gr_no_overlap,
                                               genes_gr_no_overlap,
                                               verbose = TRUE)
    expect_s3_class(region_map_no_overlap, "data.frame")
    expect_equal(nrow(region_map_no_overlap), 0)
})
