regionalpcs
================

# Table of Contents

1.  [Introduction](#introduction)
2.  [Installation](#installation)
3.  [regionalpcs R Package Tutorial](#regionalpcs-r-package-tutorial)
    - 3.1 [Loading Required Packages](#loading-required-packages)
    - 3.2 [Load the Dataset](#load-the-dataset)
      - 3.2.1 [Overview](#overview)
      - 3.2.2 [Inspecting the Data](#inspecting-the-data)
    - 3.3 [Obtaining Methylation Array Probe
      Positions](#obtaining-methylation-array-probe-positions)
      - 3.3.1 [Introduction](#introduction-1)
      - 3.3.2 [Extract Probe Names and
        Positions](#extract-probe-names-and-positions)
      - 3.3.3 [Load Illumina 450k Array Probe
        Positions](#load-illuminaprobe-positions)
      - 3.3.4 [Merge Data Frames](#merge-data-frames)
      - 3.3.5 [Addressing Genome Build
        Discrepancies](#addressing-genome-build-discrepancies)
    - 3.4 [Processing and Filtering Methylation
      Data](#processing-and-filtering-methylation-data)
      - 3.4.1 [Introduction](#introduction-2)
      - 3.4.2 [Remove Low Variance CpGs](#remove-low-variance-cpgs)
      - 3.4.3 [Normalize Methylation
        Values](#normalize-methylation-values)
    - 3.5 [Summarizing Gene Region
      Types](#summarizing-gene-region-types)
      - 3.5.1 [Introduction](#introduction-3)
      - 3.5.2 [Load Gene Region
        Annotations](#load-gene-region-annotations)
      - 3.5.3 [Create a Region Map](#create-a-region-map)
        - 3.5.3.1 [Extract CpG Positions](#extract-cpg-positions)
        - 3.5.3.2 [Convert to GenomicRanges and Find
          Overlaps](#convert-to-genomicranges-and-find-overlaps)
      - 3.5.4 [Summarizing Gene Regions with Regional Principal
        Components](#summarizing-gene-regions)
        - 3.5.4.1 [Compute Regional PCs](#compute-regional-pcs)
        - 3.5.4.2 [Inspecting the Output](#inspecting-the-output)
        - 3.5.4.3 [Extracting and Viewing Regional
          PCs](#extracting-and-viewing-regional-pcs)
        - 3.5.4.4 [Understanding the
          Results](#understanding-the-results)

# regionalpcs

Tiffany Eulalio

The `regionalpcs` package aims to address the challenge of summarizing
and interpreting DNA methylation data at a regional level. Traditional
methods of analysis may not capture the biological complexity of
methylation patterns, potentially leading to less accurate or less
meaningful interpretations. This package introduces the concept of
regional principal components (rPCs) as a tool for capturing more
biologically relevant signals in DNA methylation data. By using rPCs,
researchers can gain new insights into complex interactions and effects
in methylation data that might otherwise be missed.

# Installation

You can install the regionalpcs package from Bioconductor using the
following command:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("regionalpcs")
```

You can install the development version of regionalpcs from GitHub with:

``` r
# install devtool package if needed
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")

# download the regionalpcs package
devtools::install_github("tyeulalio/regionalpcs")
```

# `regionalpcs` R Package Tutorial

## Loading Required Packages

This tutorial depends on several Bioconductor packages. These packages
should be loaded at the beginning of the analysis.

``` r
library(regionalpcs)
library(RNOmni)
library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(liftOver)
library(magrittr)
library(tidyr)
library(tibble)
library(dplyr)
```

Here, we load the regionalpcs package, which is the main package we’ll
be using in this tutorial. We also load RNOmni, which provides
normalization functions, GenomicRanges, which provides tools for working
with genomic intervals, and tidyverse, which provides a suite of tools
for data manipulation and visualization.

It’s important to note that you need to have these packages installed on
your machine before loading them. You can install them using the
install.packages() function in R.

Once the packages are loaded, we can start using the functions provided
by each package.

## Load the dataset

### Overview

The betas dataset in the `regionalpcs` package is a subset of 450k array
methylation data from TCGA, containing 293 methylation sites and 300
samples. We’ll load this data into our R session and explore its
structure.

``` r
data("betas", package = "regionalpcs")
```

### Inspecting the Data

We can take a quick look at the dimensions of the dataset and the first
few rows to understand its structure.

``` r
head(betas)[, 1:3]
#>                                     TCGA-EJ-7781-11A TCGA-BH-A1FE-11B
#> chr16_53434200_53434201_cg00000029        0.20361112       0.13654490
#> chr15_22838620_22838621_cg00000622        0.01311223       0.01024075
#> chr1_166989202_166989203_cg00001349       0.72180841       0.74037266
#> chr8_119416178_119416179_cg00002464       0.05881476       0.05834758
#> chr6_169751536_169751537_cg00005543       0.01868565       0.01808436
#> chr12_52069532_52069533_cg00006122        0.05748535       0.06136361
#>                                     TCGA-BH-A0C3-11A
#> chr16_53434200_53434201_cg00000029        0.12996001
#> chr15_22838620_22838621_cg00000622        0.01847991
#> chr1_166989202_166989203_cg00001349       0.76361838
#> chr8_119416178_119416179_cg00002464       0.06189946
#> chr6_169751536_169751537_cg00005543       0.02085631
#> chr12_52069532_52069533_cg00006122        0.05805733
dim(betas)
#> [1] 293 300
```

Note that the row names are CpG IDs and genomic positions, and the
columns contain methylation beta values ranging from 0 to 1 for
individual samples.

## Obtaining Methylation Array Probe Positions

### Introduction

To perform accurate and informative analyses on methylation array data,
it is critical to have precise genomic positions for each probe. The
`IlluminaHumanMethylation450kanno.ilmn12.hg19` package contains
annotations for 450k methylation arrays, which can be utilized for this
purpose. This section will walk you through the steps to associate each
probe in your dataset with its genomic position.

### Extract Probe Names and Positions

First, we’ll extract the probe names from the betas data frame and use
regular expressions to separate the CpG identifier from its genomic
position.

``` r
# Extract probe names and CpG positions from row names of 'betas'
cpg_df <- data.frame(cpgs = rownames(betas)) %>%
    separate(cpgs,
        into = c("cpg_pos", "probe"), sep = "_(?=[^_]+$)",
        extra = "merge"
    )
head(cpg_df)
#>                    cpg_pos      probe
#> 1  chr16_53434200_53434201 cg00000029
#> 2  chr15_22838620_22838621 cg00000622
#> 3 chr1_166989202_166989203 cg00001349
#> 4 chr8_119416178_119416179 cg00002464
#> 5 chr6_169751536_169751537 cg00005543
#> 6  chr12_52069532_52069533 cg00006122
```

### Load Illumina 450k Array Probe Positions

Next, let’s load the Illumina 450k array probe positions for further
annotation.

``` r
data(Locations)
probe_positions <- data.frame(Locations)
head(probe_positions)
#>             chr      pos strand
#> cg00050873 chrY  9363356      -
#> cg00212031 chrY 21239348      -
#> cg00213748 chrY  8148233      -
#> cg00214611 chrY 15815688      -
#> cg00455876 chrY  9385539      -
#> cg01707559 chrY  6778695      +
```

### Merge Data Frames

Now, we merge the extracted probe names with the Illumina 450k array
probe positions.

``` r
formatted_probe_positions <- probe_positions %>%
    rownames_to_column("probe")

new_cpg_df <- cpg_df %>%
    left_join(formatted_probe_positions, by = "probe")
head(new_cpg_df)
#>                    cpg_pos      probe   chr       pos strand
#> 1  chr16_53434200_53434201 cg00000029 chr16  53468112      +
#> 2  chr15_22838620_22838621 cg00000622 chr15  23034447      +
#> 3 chr1_166989202_166989203 cg00001349  chr1 166958439      -
#> 4 chr8_119416178_119416179 cg00002464  chr8 120428418      +
#> 5 chr6_169751536_169751537 cg00005543  chr6 170151632      +
#> 6  chr12_52069532_52069533 cg00006122 chr12  52463316      +
```

### Addressing Genome Build Discrepancies

It’s critical to ensure that the genome builds match across datasets. In
this example, we’ll use the `GenomicRanges` and `liftOver` packages to
convert the genomic positions from hg19 to hg38. Here’s a quick example
on how to lift over positions from one build to another. **Always ensure
that you are working with the correct genome build and that the build
matches across all your datasets, or else you will run into big
issues!**

We need a chain file to lift the genomic positions. The chain file is an
annotation file that links the positions between the genome builds. You
can download this file from the (UCSC golden path download
site)\[<https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/>\]. Be
sure to download the file that maps between the appropriate builds.
We’ll be mapping from hg19 to hg38. We’ve included the chain used in
this analysis as a part of the regionalpcs package, which can be
accessed in the “extdata” folder as shown in the code below.

``` r
# Convert hg19 positions into a GenomicRanges object
hg19_pos <- new_cpg_df %>%
    select("chr", "pos", "strand", "probe") %>%
    mutate(start = pos, end = pos + 1)

hg19_pos_gr <- makeGRangesFromDataFrame(hg19_pos, keep.extra.columns = TRUE)

# Load chain file and liftOver positions
chain_file <- system.file("extdata", "hg19ToHg38.over.chain",
    package = "regionalpcs"
)
print(paste("Using chain file for liftOver", chain_file))
#> [1] "Using chain file for liftOver C:/Users/tyeul/AppData/Local/R/win-library/4.3/regionalpcs/extdata/hg19ToHg38.over.chain"
print(file.exists(chain_file))
#> [1] TRUE
chain <- import.chain(chain_file)

hg38_pos <- liftOver(hg19_pos_gr, chain) %>%
    as.data.frame()

# Merge the lifted positions back to the original data frame
formatted_hg38 <- hg38_pos %>%
    select(chrom_hg38 = seqnames, pos_hg38 = start, probe)

lifted_cpg_df <- new_cpg_df %>%
    left_join(formatted_hg38, by = "probe")
head(lifted_cpg_df)
#>                    cpg_pos      probe   chr       pos strand chrom_hg38
#> 1  chr16_53434200_53434201 cg00000029 chr16  53468112      +      chr16
#> 2  chr15_22838620_22838621 cg00000622 chr15  23034447      +      chr15
#> 3 chr1_166989202_166989203 cg00001349  chr1 166958439      -       chr1
#> 4 chr8_119416178_119416179 cg00002464  chr8 120428418      +       chr8
#> 5 chr6_169751536_169751537 cg00005543  chr6 170151632      +       chr6
#> 6  chr12_52069532_52069533 cg00006122 chr12  52463316      +      chr12
#>    pos_hg38
#> 1  53434200
#> 2  22838620
#> 3 166989202
#> 4 119416178
#> 5 169751536
#> 6  52069532
```

Now that we have accurate genomic positions for each probe and have
harmonized genome builds, we can proceed with preprocessing the
methylation data.

## Processing and Filtering Methylation Data

### Introduction

Before conducting downstream analyses, it is essential to preprocess and
clean the methylation data. In this section, we’ll walk you through the
steps to remove low variance CpGs and normalize the methylation beta
values.

### Remove Low Variance CpGs

Firstly, we aim to filter out low variance CpGs. Variability is a
crucial factor, as low variance CpGs may not provide much information
for downstream analyses.

In this section, we’ll remove low variance CpGs and normalize our
methylation beta values using the inverse normal transform.

``` r
# Remove CpGs with zero variance
var_betas <- betas[apply(betas, 1, var, na.rm = TRUE) != 0, ] %>%
    na.omit()
dim(var_betas)
#> [1] 293 300
```

We only remove CpGs that have zero variance in this example. You can
adjust this threshold according to the requirements of your specific
analysis.

### Normalize Methylation Values

Methylation data often exhibit heteroscedasticity. Therefore, we’ll
normalize the beta values using inverse normal transformation. For this,
we’ll use the `RankNorm` function from the `RNOmni` package.

``` r
# Apply inverse normal transformation to methylation beta values
int_meth <- apply(var_betas, 1, RankNorm) %>%
    t() %>%
    as.data.frame()
```

After these preprocessing steps, you will have a dataset ready for
downstream analysis with the `regionalpcs` package. We’ll cover how to
perform these analyses in subsequent sections of this tutorial.

## Summarizing Gene Region Types

### Introduction

Gene regions are significant functional units of the genome, such as
promoters, gene bodies, and intergenic regions. We’ll focus on
summarizing these regions to prepare for downstream analyses. We will
use the `regionalpcs` package to perform these tasks.

### Load Gene Region Annotations

First, let’s load the gene region annotations. Make sure to align the
genomic builds of your annotations and methylation data.

**All annotations included with the `regionalpcs` package are in build
hg38.**

``` r
# Load the gene region annotation file
data("gene_annots", package = "regionalpcs")
head(gene_annots)
#> # A tibble: 6 × 16
#>   seqnames   start     end width strand tx_id             type   gencode_gene_id
#>   <chr>      <dbl>   <dbl> <dbl> <chr>  <chr>             <chr>  <chr>          
#> 1 chr1      922928  923927  1000 +      ENST00000420190.6 hg38_… ENSG0000018763…
#> 2 chr1      959584  960583  1000 +      ENST00000338591.8 hg38_… ENSG0000018796…
#> 3 chr1      965482  966481  1000 +      ENST00000379410.8 hg38_… ENSG0000018758…
#> 4 chr1     1000138 1001137  1000 +      ENST00000624697.4 hg38_… ENSG0000018760…
#> 5 chr1     1019120 1020119  1000 +      ENST00000379370.7 hg38_… ENSG0000018815…
#> 6 chr1     1172903 1173902  1000 +      ENST00000379290.5 hg38_… ENSG0000016257…
#> # ℹ 8 more variables: gencode_gene_type <chr>, gencode_gene_name <chr>,
#> #   transcript_type <chr>, transcript_name <chr>,
#> #   transcript_support_level <dbl>, tag <chr>, is_canonical <lgl>,
#> #   gencode_region <chr>
```

The `gene_annots` dataset includes annotations for various gene regions.

### Create a Region Map

Before summarizing gene regions using `compute_regional_pcs`, we need to
create a region map that assigns CpGs to gene regions. This map enables
us to identify which CpGs fall into each gene region.

#### Extract CpG Positions

Start by extracting the CpG positions from your methylation data frame’s
row names.

``` r
head(int_meth)[1:4]
#>                                     TCGA-EJ-7781-11A TCGA-BH-A1FE-11B
#> chr16_53434200_53434201_cg00000029       -0.38948253       -0.9465265
#> chr15_22838620_22838621_cg00000622       -0.03757696       -1.1171114
#> chr1_166989202_166989203_cg00001349      -0.87083034       -0.6274087
#> chr8_119416178_119416179_cg00002464       0.13818831        0.1045460
#> chr6_169751536_169751537_cg00005543       0.38049184        0.2488238
#> chr12_52069532_52069533_cg00006122        0.27474294        0.6274087
#>                                     TCGA-BH-A0C3-11A TCGA-E9-A1N6-11A
#> chr16_53434200_53434201_cg00000029        -1.1171114       -1.4370361
#> chr15_22838620_22838621_cg00000622         1.2160129       -0.6376059
#> chr1_166989202_166989203_cg00001349       -0.3894825        0.3894825
#> chr8_119416178_119416179_cg00002464        0.2402219       -0.2574441
#> chr6_169751536_169751537_cg00005543        0.9080280       -1.5657713
#> chr12_52069532_52069533_cg00006122         0.3271598       -0.1129440
# Extract CpG information
cpg_info <- data.frame(cpg_id = rownames(int_meth)) %>%
    separate(cpg_id,
        into = c("chrom", "start", "end", "cpg_name"),
        sep = "_", remove = FALSE
    )
head(cpg_info)
#>                                cpg_id chrom     start       end   cpg_name
#> 1  chr16_53434200_53434201_cg00000029 chr16  53434200  53434201 cg00000029
#> 2  chr15_22838620_22838621_cg00000622 chr15  22838620  22838621 cg00000622
#> 3 chr1_166989202_166989203_cg00001349  chr1 166989202 166989203 cg00001349
#> 4 chr8_119416178_119416179_cg00002464  chr8 119416178 119416179 cg00002464
#> 5 chr6_169751536_169751537_cg00005543  chr6 169751536 169751537 cg00005543
#> 6  chr12_52069532_52069533_cg00006122 chr12  52069532  52069533 cg00006122
```

#### Convert to GenomicRanges and Find Overlaps

Next, we’ll use the `GenomicRanges` package to find overlaps between
CpGs and gene regions.

``` r
# Convert to GenomicRanges
cpg_gr <- makeGRangesFromDataFrame(cpg_info, keep.extra.columns = TRUE)
annots_gr <- makeGRangesFromDataFrame(gene_annots, keep.extra.columns = TRUE)

# Find overlaps between the two GRanges objects
overlaps <- findOverlaps(query = cpg_gr, subject = annots_gr) %>%
    as.data.frame()
head(overlaps)
#>   queryHits subjectHits
#> 1         1       12678
#> 2         2       11904
#> 3         3         679
#> 4         4        7360
#> 5         5        5877
#> 6         6       10306

# Match overlaps
matched_cpg <- cpg_gr[overlaps$queryHits, ] %>%
    as.data.frame() %>%
    select(cpg_id)

# Select overlapped rows and just keep the columns we need
matched_annots <- annots_gr[overlaps$subjectHits, ] %>%
    as.data.frame() %>%
    select(gencode_gene_id)

# Combine the matched CpGs and gene annotations to form the region map
region_map <- cbind(matched_annots, matched_cpg)
head(region_map)
#>      gencode_gene_id                              cpg_id
#> 1 ENSG00000103479.16  chr16_53434200_53434201_cg00000029
#> 2 ENSG00000140157.14  chr15_22838620_22838621_cg00000622
#> 3 ENSG00000143194.13 chr1_166989202_166989203_cg00001349
#> 4  ENSG00000136999.5 chr8_119416178_119416179_cg00002464
#> 5 ENSG00000130023.16 chr6_169751536_169751537_cg00005543
#> 6 ENSG00000123395.14  chr12_52069532_52069533_cg00006122
length(unique(region_map$gencode_gene_id))
#> [1] 52
```

With these steps, you’ll have a region map that assigns CpGs to specific
gene regions, which can be essential for downstream analyses.

<a name="summarizing-gene-regions"></a>

### Summarizing Gene Regions with Regional Principal Components

In this final section, we’ll summarize gene regions using Principal
Components (PCs) to capture the maximum variation. We’ll utilize the
`compute_regional_pcs` function from the `regionalpcs` package for this.

#### Compute Regional PCs

Let’s calculate the regional PCs using a subset of our gene regions for
demonstration purposes.

``` r
# Display head of region map
head(region_map)
#>      gencode_gene_id                              cpg_id
#> 1 ENSG00000103479.16  chr16_53434200_53434201_cg00000029
#> 2 ENSG00000140157.14  chr15_22838620_22838621_cg00000622
#> 3 ENSG00000143194.13 chr1_166989202_166989203_cg00001349
#> 4  ENSG00000136999.5 chr8_119416178_119416179_cg00002464
#> 5 ENSG00000130023.16 chr6_169751536_169751537_cg00005543
#> 6 ENSG00000123395.14  chr12_52069532_52069533_cg00006122

# Subset the region map
sub_region_map <- region_map %>%
    filter(gencode_gene_id %in% unique(region_map$gencode_gene_id)[1:1000])

# Compute regional PCs
res <- compute_regional_pcs(int_meth, sub_region_map)
#> Using Gavish-Donoho method
```

#### Inspecting the Output

The function returns a list containing multiple elements. Let’s first
look at what these elements are.

``` r
# Inspect the output list elements
names(res)
#> [1] "regional_pcs"     "percent_variance" "loadings"         "info"
```

#### Extracting and Viewing Regional PCs

The first element (`res$regional_pcs`) is a data frame containing the
calculated regional PCs.

``` r
# Extract regional PCs
regional_pcs <- res$regional_pcs
head(regional_pcs)[1:4]
#>                        TCGA-EJ-7781-11A TCGA-BH-A1FE-11B TCGA-BH-A0C3-11A
#> ENSG00000103479.16-PC1       -0.8210157      -1.60659103       -1.0541882
#> ENSG00000140157.14-PC1       -0.1392372      -1.47241611        1.3317941
#> ENSG00000143194.13-PC1       -1.7895014      -3.07334567       -3.1270943
#> ENSG00000136999.5-PC1         0.4498687       0.11381324       -0.3555689
#> ENSG00000130023.16-PC1       -0.1838917       0.05055839       -0.8546051
#> ENSG00000123395.14-PC1       -1.1961240       0.36432710        0.7234163
#>                        TCGA-E9-A1N6-11A
#> ENSG00000103479.16-PC1       -1.6169547
#> ENSG00000140157.14-PC1       -0.6519472
#> ENSG00000143194.13-PC1       -1.0122026
#> ENSG00000136999.5-PC1         1.2443734
#> ENSG00000130023.16-PC1        2.1260186
#> ENSG00000123395.14-PC1        1.3592205
```

#### Understanding the Results

The output is a data frame with regional PCs for each region as rows and
our samples as columns. This is our new representation of methylation
values, now on a gene regional PC scale. We can feed these into
downstream analyses as is.

The number of regional PCs representing each gene region was determined
by the Gavish-Donoho method. This method allows us to identify PCs that
capture actual signal in our data and not the noise that is inherent in
any dataset. To explore alternative methods, we can change the
`pc_method` parameter.

``` r
# Count the number of unique gene regions and PCs
regions <- data.frame(gene_pc = rownames(regional_pcs)) %>%
    separate(gene_pc, into = c("gene", "pc"), sep = "-")
head(regions)
#>                 gene  pc
#> 1 ENSG00000103479.16 PC1
#> 2 ENSG00000140157.14 PC1
#> 3 ENSG00000143194.13 PC1
#> 4  ENSG00000136999.5 PC1
#> 5 ENSG00000130023.16 PC1
#> 6 ENSG00000123395.14 PC1

# number of genes that have been summarized
length(unique(regions$gene))
#> [1] 52

# how many of each PC did we get
table(regions$pc)
#> 
#> PC1 
#>  52
```

We have summarized each of our genes using just one PC. The number of
PCs depends on three main factors: the number of samples, the number of
CpGs in the gene region, and the noise in the methylation data.

By default, the `compute_regional_pcs` function uses the Gavish-Donoho
method. However, we can also use the Marcenko-Pasteur method by setting
the `pc_method` parameter:

``` r
# Using Marcenko-Pasteur method
mp_res <- compute_regional_pcs(int_meth, sub_region_map, pc_method = "mp")
#> Using Marchenko-Pastur method

# select the regional pcs
mp_regional_pcs <- mp_res$regional_pcs

# separate the genes from the pc numbers
mp_regions <- data.frame(gene_pc = rownames(mp_regional_pcs)) %>%
    separate(gene_pc, into = c("gene", "pc"), sep = "-")
head(mp_regions)
#>                 gene  pc
#> 1 ENSG00000103479.16 PC1
#> 2 ENSG00000103479.16 PC2
#> 3 ENSG00000140157.14 PC1
#> 4 ENSG00000140157.14 PC2
#> 5 ENSG00000143194.13 PC1
#> 6 ENSG00000143194.13 PC2

# number of genes that have been summarized
length(unique(mp_regions$gene))
#> [1] 52

# how many of each PC did we get
table(mp_regions$pc)
#> 
#> PC1 PC2 
#>  52  26
```

The Marcenko-Pasteur and the Gavish-Donoho methods are both based on
Random Matrix Theory, and they aim to identify the number of significant
PCs that capture the true signal in the data and not just the noise.
However, these methods differ in how they select the number of
significant PCs. The Marcenko-Pasteur method typically selects more PCs
to represent a gene region compared to the Gavish-Donoho method. This
may be due to the different ways in which the two methods estimate the
noise level in the data.

Ultimately, the choice between the two methods depends on the specific
needs and goals of the analysis. The Gavish-Donoho method tends to
provide more conservative results, while the Marcenko-Pasteur method may
capture more of the underlying signal in the data. Researchers should
carefully consider their objectives and the characteristics of their
data when selecting a method.

# Session Information

``` r
sessionInfo()
#> R version 4.3.1 (2023-06-16 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19045)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/Los_Angeles
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] dplyr_1.1.2                                       
#>  [2] tibble_3.2.1                                      
#>  [3] tidyr_1.3.0                                       
#>  [4] magrittr_2.0.3                                    
#>  [5] liftOver_1.25.0                                   
#>  [6] Homo.sapiens_1.3.1                                
#>  [7] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2           
#>  [8] org.Hs.eg.db_3.17.0                               
#>  [9] GO.db_3.17.0                                      
#> [10] OrganismDbi_1.43.0                                
#> [11] GenomicFeatures_1.53.1                            
#> [12] AnnotationDbi_1.63.2                              
#> [13] rtracklayer_1.61.1                                
#> [14] gwascat_2.33.0                                    
#> [15] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1
#> [16] minfi_1.47.0                                      
#> [17] bumphunter_1.43.0                                 
#> [18] locfit_1.5-9.8                                    
#> [19] iterators_1.0.14                                  
#> [20] foreach_1.5.2                                     
#> [21] Biostrings_2.69.2                                 
#> [22] XVector_0.41.1                                    
#> [23] SummarizedExperiment_1.31.1                       
#> [24] Biobase_2.61.0                                    
#> [25] MatrixGenerics_1.13.1                             
#> [26] matrixStats_1.0.0                                 
#> [27] GenomicRanges_1.53.1                              
#> [28] GenomeInfoDb_1.37.2                               
#> [29] IRanges_2.35.2                                    
#> [30] S4Vectors_0.39.1                                  
#> [31] BiocGenerics_0.47.0                               
#> [32] RNOmni_1.0.1                                      
#> [33] regionalpcs_0.99.1                                
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.3.1             BiocIO_1.11.0            
#>   [3] bitops_1.0-7              filelock_1.0.2           
#>   [5] preprocessCore_1.63.1     graph_1.79.0             
#>   [7] XML_3.99-0.14             lifecycle_1.0.3          
#>   [9] lattice_0.21-8            MASS_7.3-60              
#>  [11] base64_2.0.1              scrime_1.3.5             
#>  [13] limma_3.57.7              rmarkdown_2.24           
#>  [15] yaml_2.3.7                doRNG_1.8.6              
#>  [17] askpass_1.1               cowplot_1.1.1            
#>  [19] DBI_1.1.3                 RColorBrewer_1.1-3       
#>  [21] abind_1.4-5               zlibbioc_1.47.0          
#>  [23] quadprog_1.5-8            purrr_1.0.2              
#>  [25] RCurl_1.98-1.12           VariantAnnotation_1.47.1 
#>  [27] rappdirs_0.3.3            GenomeInfoDbData_1.2.10  
#>  [29] RMTstat_0.3.1             ggrepel_0.9.3            
#>  [31] irlba_2.3.5.1             genefilter_1.83.1        
#>  [33] dqrng_0.3.0               annotate_1.79.0          
#>  [35] DelayedMatrixStats_1.23.4 codetools_0.2-19         
#>  [37] DelayedArray_0.27.10      xml2_1.3.5               
#>  [39] tidyselect_1.2.0          ScaledMatrix_1.9.1       
#>  [41] beanplot_1.3.1            BiocFileCache_2.9.1      
#>  [43] illuminaio_0.43.0         GenomicAlignments_1.37.0 
#>  [45] multtest_2.57.0           survival_3.5-7           
#>  [47] tools_4.3.1               progress_1.2.2           
#>  [49] Rcpp_1.0.11               glue_1.6.2               
#>  [51] SparseArray_1.1.11        xfun_0.40                
#>  [53] HDF5Array_1.29.3          withr_2.5.0              
#>  [55] BiocManager_1.30.22       fastmap_1.1.1            
#>  [57] rhdf5filters_1.13.5       fansi_1.0.4              
#>  [59] openssl_2.1.0             rsvd_1.0.5               
#>  [61] digest_0.6.33             R6_2.5.1                 
#>  [63] colorspace_2.1-0          biomaRt_2.57.1           
#>  [65] RSQLite_2.3.1             utf8_1.2.3               
#>  [67] generics_0.1.3            data.table_1.14.8        
#>  [69] prettyunits_1.1.1         httr_1.4.7               
#>  [71] S4Arrays_1.1.5            pkgconfig_2.0.3          
#>  [73] gtable_0.3.3              blob_1.2.4               
#>  [75] siggenes_1.75.0           htmltools_0.5.6          
#>  [77] RBGL_1.77.1               scales_1.2.1             
#>  [79] png_0.1-8                 knitr_1.43               
#>  [81] rstudioapi_0.15.0         reshape2_1.4.4           
#>  [83] tzdb_0.4.0                rjson_0.2.21             
#>  [85] nlme_3.1-163              curl_5.0.2               
#>  [87] cachem_1.0.8              rhdf5_2.45.1             
#>  [89] stringr_1.5.0             restfulr_0.0.15          
#>  [91] GEOquery_2.69.0           pillar_1.9.0             
#>  [93] grid_4.3.1                reshape_0.8.9            
#>  [95] vctrs_0.6.3               BiocSingular_1.17.1      
#>  [97] beachmat_2.17.15          dbplyr_2.3.3             
#>  [99] xtable_1.8-4              evaluate_0.21            
#> [101] readr_2.1.4               cli_3.6.1                
#> [103] compiler_4.3.1            Rsamtools_2.17.0         
#> [105] rlang_1.1.1               crayon_1.5.2             
#> [107] rngtools_1.5.2            nor1mix_1.3-0            
#> [109] mclust_6.0.0              plyr_1.8.8               
#> [111] stringi_1.7.12            BiocParallel_1.35.3      
#> [113] munsell_0.5.0             PCAtools_2.13.0          
#> [115] Matrix_1.6-1              BSgenome_1.69.0          
#> [117] hms_1.1.3                 sparseMatrixStats_1.13.4 
#> [119] bit64_4.0.5               ggplot2_3.4.3            
#> [121] Rhdf5lib_1.23.0           KEGGREST_1.41.0          
#> [123] statmod_1.5.0             memoise_2.0.1            
#> [125] snpStats_1.51.0           bit_4.0.5
```
