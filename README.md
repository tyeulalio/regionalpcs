regionalpcs
================

# Table of Contents

1.  [Introduction](#introduction)
2.  [Repository Contents](#repository-contents)
3.  [System Requirements](#system-requirements)
4.  [Installation Guide](#installation)
5.  [Demo](#demo)

# Introduction

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

# Repository Contents

- `R/`: Contains the source code for the project, written in R. This
  directory includes all the scripts and functions that constitute the
  core functionality of the package.

- `inst/`: Stores files that are retained post-installation of the R
  package. This includes additional data, documentation, or scripts that
  users might find useful when working with the package.

- `man/`: Contains the manual pages for the package. These documentation
  files are accessible within an R session using the `help()` function,
  providing users with detailed information on package functions and
  usage.

- `tests/`: Houses R unit tests developed with the `testthat` package.
  These tests are designed to automatically verify that the package
  functions correctly under various conditions, ensuring reliability and
  stability.

- `vignettes/`: Provides R Markdown vignettes that offer comprehensive
  tutorials and examples on how to use the package. These vignettes are
  converted into HTML help pages accessible from within an R session,
  serving as a valuable resource for users to learn about the package
  features and functionalities.

# System Requirements

## Hardware Requirements

The `regionalpcs` package is designed to function efficiently on a
standard computer setup. The specific RAM requirement depends on the
scale of the analysis defined by the user. Below are our recommendations
for minimal and optimal performance configurations:

- **Minimal Configuration**: A computer with at least 2 GB of RAM is
  required for basic operation.
- **Recommended Configuration**: For optimal performance, especially for
  more demanding analyses, we recommend the following specifications:
  - **RAM**: 16 GB or more
  - **CPU**: 4 or more cores, with a clock speed of 3.3 GHz per core or
    faster

**Runtime Benchmarks**: The reported runtimes are based on tests
conducted on a system equipped with 64 GB RAM, an 8-core CPU @ 3.60 GHz,
and an internet connection speed of 229 Mbps.

## Software Requirements

### Operating System Compatibility

While the development version of the `regionalpcs` package is primarily
tested on Windows platforms, we aim for broad compatibility across major
operating systems.

Our Bioconductor package
[`regionalpcs`](https://bioconductor.org/packages/release/bioc/html/regionalpcs.html)
has been tested with Windows, Mac, and Linux operating systems. Here are
the details regarding the tested systems:

- **Linux**: Ubuntu 22.04.03 LTS / x86_64
- **Mac OSX**: macOS 12.7.1 Monterey / x86_64
- **Windows**: 10 Pro / x64

### R and Package Dependencies

To install and run the `regionalpcs` package, the following software
requirements must be met:

- **R Version**: The package requires R version 4.3.0 or higher. Ensure
  that your R installation is up to date before proceeding with the
  installation of `regionalpcs`.
- **Dependencies**: Additional R packages from CRAN and possibly
  Bioconductor are required. Users will be prompted to install any
  missing dependencies during the package installation process. Required
  packages are:

``` r
dplyr
PCAtools 
tibble
GenomicRanges
```

Please refer to the package documentation for a detailed list of
dependencies and instructions for setting up the required software
environment.

# Installation Guide

You can install the regionalpcs package from Bioconductor using the
following command:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("regionalpcs")
```

which will install in about 30 seconds on a machine with the recommended
specs.

You can install the development version of regionalpcs from GitHub with:

``` r
# install devtool package if needed
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")

# download the regionalpcs package
devtools::install_github("tyeulalio/regionalpcs")
```

# Demonstration

Explore the functionalities of the `regionalpcs` package with our
interactive tutorials provided as vignettes. These vignettes offer
step-by-step guidance on using the package’s main features and are
designed to help you get started quickly.

### Accessing the Vignettes

To start the tutorials, ensure that the `regionalpcs` package is
installed and loaded into your R session. You can then access the
vignettes directly in R with the following commands:

``` r
# Load the regionalpcs package
library(regionalpcs)

# Open the main vignette
vignette('regionalpcs-introduction')
```

### Online Access

Alternatively, for access to a browser-friendly version, visit the
[`regionalpcs` Bioconductor
page](https://bioconductor.org/packages/release/bioc/html/regionalpcs.html).
Here, you’ll find the vignettes available in HTML and R formats.

### Tutorial Duration

The primary vignette is concise and informative, designed to provide a
comprehensive overview within approximately 20 seconds. This makes it an
efficient way to familiarize yourself with the package’s capabilities
and start applying them to your data analysis projects.

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
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.3.1    fastmap_1.1.1     cli_3.6.1         tools_4.3.1      
#>  [5] htmltools_0.5.6   rstudioapi_0.16.0 yaml_2.3.7        rmarkdown_2.26   
#>  [9] knitr_1.46        xfun_0.43         digest_0.6.33     rlang_1.1.1      
#> [13] evaluate_0.23
```
