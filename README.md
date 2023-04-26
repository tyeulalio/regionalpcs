
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regionalpcs

<!-- badges: start -->
<!-- badges: end -->

The goal of regionalpcs is to provide functions to easily summarize DNA
methylation data on a regional level using regional principal components
(rPCs). The rPCs help to capture more biologically-relevant signal for
downstream analyses.

## Installation

You can install the development version of regionalpcs from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tyeulalio/regionalpcs")
```

# `regionalpcs` R Package Tutorial

## Loading packages

To begin working with the *regionalpcs* R package, we need to first load
it into our R session along with some other necessary packages. We can
do this using the **`library()`** function in R. Here’s an example code
snippet to load the required packages:

``` r
library(regionalpcs)
library(RNOmni)
library(GenomicRanges)
library(tidyverse)
```

Here, we load the *regionalpcs* package, which is the main package we’ll
be using in this tutorial. We also load *RNOmni*, which provides
normalization functions, *GenomicRanges*, which provides tools for
working with genomic intervals, and *tidyverse*, which provides a suite
of tools for data manipulation and visualization.

It’s important to note that you need to have these packages installed on
your machine before loading them. You can install them using the
**`install.packages()`** function in R.

Once the packages are loaded, we can start using the functions provided
by each package.

## Loading in a dataset

In this section, we’ll load a methylation dataset collected from TCGA
using 450k methylation arrays. Before we can start any analysis, we need
to load in the data and take a look at what we’re working with.

The **`betas`** dataset included with the `regionalpcs` package contains
methylation values for 1,000 methylation sites and 300 individuals.
We’ve subset the original dataset to save space and compute time for
this vignette. To work well for downstream statistical analyses, we’ll
need to normalize the beta values by removing low variance CpGs and
applying the inverse normal transform.

**remove NA values from the actual dataset - this dataset needs to be
preprocessed**

``` r
head(betas)[1:3]
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

**add something about getting positions from cpg id’s**

We can see that the row names contain CpG IDs and positions, while the
columns contain methylation beta values, which range from 0 to 1.

Before we can start analyzing the data, we need to preprocess it. This
includes removing any missing values and normalizing the beta values.
We’ll cover these steps in more detail later in the tutorial.

## Process and filter the methylation data.

In this section, we’ll remove low variance CpGs and normalize our
methylation beta values using the inverse normal transform.

``` r
# remove low variance cpgs
var_betas <- betas[apply(betas, 1, var, na.rm=TRUE) != 0,] %>%
  na.omit()
dim(var_betas)
#> [1] 293 300
```

We’ll be lenient here and only remove CpGs that have zero variance.
However, we can come back and change this threshold later if needed.

Now, let’s normalize our methylation values to deal with the
heteroscedasticity present in methylation data. We’ll apply the inverse
normal transform using functions from the **`RNOmni`** package by
applying the **`RankNorm`** function to each row of our data frame.

``` r
# inverse normal transform our methylation beta values
int_meth <- apply(var_betas, 1, RankNorm) %>%
  t() %>%
  as.data.frame()
```

With our filtered and normalized data, we can now summarize region types
using the **`regionalpcs`** package.

## Summarizing gene region types

### Loading gene region annotations

We’ll start by loading gene region annotations. Make sure that these
annotations are using the same genomic reference build (GrCh37, GrCh38)
as your CpG annotations.

``` r
#  load gene region annotation file
head(gene_annots)
#> # A tibble: 6 x 16
#>   seqna~1  start    end width strand tx_id type  genco~2 genco~3 genco~4 trans~5
#>   <chr>    <dbl>  <dbl> <dbl> <chr>  <chr> <chr> <chr>   <chr>   <chr>   <chr>  
#> 1 chr1    9.23e5 9.24e5  1000 +      ENST~ hg38~ ENSG00~ protei~ SAMD11  protei~
#> 2 chr1    9.60e5 9.61e5  1000 +      ENST~ hg38~ ENSG00~ protei~ KLHL17  protei~
#> 3 chr1    9.65e5 9.66e5  1000 +      ENST~ hg38~ ENSG00~ protei~ PLEKHN1 protei~
#> 4 chr1    1.00e6 1.00e6  1000 +      ENST~ hg38~ ENSG00~ protei~ ISG15   protei~
#> 5 chr1    1.02e6 1.02e6  1000 +      ENST~ hg38~ ENSG00~ protei~ AGRN    protei~
#> 6 chr1    1.17e6 1.17e6  1000 +      ENST~ hg38~ ENSG00~ protei~ TTLL10  protei~
#> # ... with 5 more variables: transcript_name <chr>,
#> #   transcript_support_level <dbl>, tag <chr>, is_canonical <lgl>,
#> #   gencode_region <chr>, and abbreviated variable names 1: seqnames,
#> #   2: gencode_gene_id, 3: gencode_gene_type, 4: gencode_gene_name,
#> #   5: transcript_type
```

The **`gene_annots`** dataset contains annotations for different gene
regions, such as promoters, gene bodies, and intergenic regions.

### Creating a region map

Before summarizing gene regions using **`compute_regional_pcs`**, we
need to create a region map that assigns CpGs to gene regions. This map
enables us to identify which CpGs fall into each gene region, and the
function **`make_region_map`** from the *regionalpcs* package automates
this process for us.

To create the region map, we need genomic positions for our CpGs, which
we can obtain from the row names of our methylation data frame.

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

# Get CpG positions from methylation data frame row names
cpg_info = data.frame(cpg_id=rownames(int_meth)) %>%
  # separate info from the cpg ID's
  separate(cpg_id, into=c('chrom', 'start', 'end', 'cpg_name'), sep='_', remove=FALSE)
head(cpg_info)
#>                                cpg_id chrom     start       end   cpg_name
#> 1  chr16_53434200_53434201_cg00000029 chr16  53434200  53434201 cg00000029
#> 2  chr15_22838620_22838621_cg00000622 chr15  22838620  22838621 cg00000622
#> 3 chr1_166989202_166989203_cg00001349  chr1 166989202 166989203 cg00001349
#> 4 chr8_119416178_119416179_cg00002464  chr8 119416178 119416179 cg00002464
#> 5 chr6_169751536_169751537_cg00005543  chr6 169751536 169751537 cg00005543
#> 6  chr12_52069532_52069533_cg00006122 chr12  52069532  52069533 cg00006122
```

Next, we convert our CpG information and gene region annotations to
`GenomicRanges` objects to find the overlaps between the two using
`findOverlaps.`

``` r
# Convert CpG info and gene annotations to GenomicRanges objects
cpg_gr <- makeGRangesFromDataFrame(cpg_info, keep.extra.columns = TRUE)
annots_gr <- makeGRangesFromDataFrame(gene_annots, keep.extra.columns = TRUE)

# Find overlaps between the two GRanges objects
overlaps <- findOverlaps(query=cpg_gr, subject=annots_gr) %>%
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
matched_cpg <- cpg_gr[overlaps$queryHits,] %>%
  as.data.frame() %>%
  select(cpg_id)

# Select overlapped rows and just keep the columns we need
matched_annots <- annots_gr[overlaps$subjectHits,] %>%
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

We can see the resulting region map, which summarizes how our CpGs are
assigned to gene regions.

### Summarize gene regions

We’ll use the `compute_regional_pcs` function to get regional PCs using
default settings.

``` r
head(region_map)
#>      gencode_gene_id                              cpg_id
#> 1 ENSG00000103479.16  chr16_53434200_53434201_cg00000029
#> 2 ENSG00000140157.14  chr15_22838620_22838621_cg00000622
#> 3 ENSG00000143194.13 chr1_166989202_166989203_cg00001349
#> 4  ENSG00000136999.5 chr8_119416178_119416179_cg00002464
#> 5 ENSG00000130023.16 chr6_169751536_169751537_cg00005543
#> 6 ENSG00000123395.14  chr12_52069532_52069533_cg00006122
sub_region_map <- region_map %>%
  filter(gencode_gene_id %in% unique(region_map$gencode_gene_id)[1:1000])

res <- compute_regional_pcs(int_meth, sub_region_map)
#> [1] "Using Gavish-Donoho method"
```

The `compute_regional_pcs` function returns a list containing important
elements from our regional PC summary method. To understand the output,
we can inspect the names of the list elements:

``` r
names(res)
#> [1] "regional_pcs"     "percent_variance" "loadings"         "info"
```

The first element of this list is the regional PCs. These are the
summaries for each region contained in the input `region_map.` To select
and inspect the regional PCs, we can do the following:

``` r
# select the regional pcs
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

The output is a data frame with regional PCs for each region as rows and
our samples as columns. This is our new representation of methylation
values, now on a gene regional PC scale. We can feed these into
downstream analyses as is.

The number of regional PCs representing each gene region was determined
by the Gavish-Donoho method. This method allows us to identify PCs that
capture actual signal in our data and not the noise that is inherent in
any dataset.

We can check the number of gene regions and regional PCs we have now:

``` r
# separate the genes from the pc numbers
regions <- data.frame(gene_pc=rownames(regional_pcs)) %>%
  separate(gene_pc, into=c('gene', 'pc'), sep='-')
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
mp_res <- compute_regional_pcs(int_meth, sub_region_map, pc_method='mp')
#> [1] "Using Marchenko-Pastur method"

# select the regional pcs
mp_regional_pcs <- mp_res$regional_pcs

# separate the genes from the pc numbers
mp_regions <- data.frame(gene_pc=rownames(mp_regional_pcs)) %>%
  separate(gene_pc, into=c('gene', 'pc'), sep='-')
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

# multi_genes <- mp_regions %>%
#   filter(pc %in% c('PC2', 'PC3', 'PC4'))
# head(multi_genes)
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
