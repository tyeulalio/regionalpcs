## just reading in methylation and gene annotation data here for now
library(readr)

# define paths here to hide them from the final document
tcga_path <- paste0("D:/montgomery_lab/deconvolution_project/making_r_package/regionalpcs_data/tcga450k/normal_exp_meth/")
gene_annots_path <- paste0("D:/montgomery_lab/deconvolution_project/making_r_package/regionalpcs_data/gene_annotations/")

# methylation data
datafile <- paste0(tcga_path, "methylation_data.rds")
betas <- readRDS(datafile)
dim(betas)
head(betas)

# keep genes that have multiple PCs using MP method
# got this list from running the tutorial on the first 1000 genes
# in region map
keep_genes <- unique(multi_genes$gene)
head(keep_genes)
length(keep_genes)
head(region_map)

# keep some of them
sub_multi <- keep_genes[1:25]

# select cpgs in the genes we're keeping
multi_map <- region_map %>%
  filter(gencode_gene_id %in% sub_multi)
length(unique(multi_map$cpg_id)) # 200 cpgs selected

# grab some other random genes
other_genes <- setdiff(unique(region_map$gencode_gene_id), sub_multi)
set.seed(71446174)
other_sample <- sample(other_genes, 20)

other_map <- region_map %>%
  filter(gencode_gene_id %in% other_sample)
head(other_map)

# how many cpgs selected all together?
length(unique(c(multi_map$cpg_id, other_map$cpg_id))) #293

# get the set of cpgs
keep_cpgs <- unique(c(multi_map$cpg_id, other_map$cpg_id))

# subset the betas
betas <- betas[keep_cpgs,1:300]
usethis::use_data(betas, overwrite = TRUE)

# annotation data
datafile <- paste0(gene_annots_path, "hg38_genes_promoters_annotations.tsv")
gene_annots <- read_tsv(datafile, show_col_types = FALSE)
head(gene_annots)

usethis::use_data(gene_annots, overwrite = TRUE)
