
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(regionalpcs)
## basic example code

head(mtcars)
#>                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
#> Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#> Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1

get_sig_pcs(t(mtcars))
#> [1] "Using Gavish-Donoho method"
#> $sig_pcs
#>                             PC1         PC2        PC3
#> Mazda RX4            -79.596425    2.132241 -2.1533361
#> Mazda RX4 Wag        -79.598570    2.147487 -2.2151237
#> Datsun 710          -133.894096   -5.057570 -2.1379503
#> Hornet 4 Drive         8.516559   44.985630  1.2337628
#> Hornet Sportabout    128.686342   30.817402  3.3434206
#> Valiant              -23.220146   35.106518 -3.2595619
#> Duster 360           159.309025  -32.259197  0.6489441
#> Merc 240D           -112.615805   39.702195 -0.4652441
#> Merc 230            -103.534591    7.513104 -1.5893097
#> Merc 280             -67.046877   -6.208536 -3.6071985
#> Merc 280C            -66.997514   -6.206387 -5.0250999
#> Merc 450SE            55.211672  -10.373509 -1.6241522
#> Merc 450SL            55.173910  -10.361893 -0.7257531
#> Merc 450SLC           55.251602  -10.370934 -2.8210156
#> Cadillac Fleetwood   242.814893   52.501758 -0.9704237
#> Lincoln Continental  236.369886   38.280788 -1.0973314
#> Chrysler Imperial    224.737944   16.111941  2.9176287
#> Fiat 128            -172.363654    6.575522  5.6488371
#> Honda Civic         -181.066911   17.783639  3.2876931
#> Toyota Corolla      -179.697852    4.188212  6.8636215
#> Toyota Corona       -121.224099   -3.345362 -3.1557542
#> Dodge Challenger      80.159386   34.983214 -1.7586589
#> AMC Javelin           67.572431   28.894067 -2.5015855
#> Camaro Z28           150.354631  -36.633575 -0.6190038
#> Pontiac Firebird     164.652522   48.239880  5.0528427
#> Fiat X1-9           -171.897231    6.643746  0.7130619
#> Porsche 914-2       -123.804988    2.033356  1.4578737
#> Lotus Europa        -137.082789  -28.675647  5.5555670
#> Ford Pantera L       159.413222  -53.318347  2.6378142
#> Ferrari Dino         -64.762396  -62.954280 -2.3829459
#> Maserati Bora        145.361703 -139.049149  1.5818981
#> Volvo 142E          -115.181783  -13.826313 -2.8335169
#> 
#> $percent_var
#>                                  PC1       PC2          PC3
#> percent_variance_explained 0.9269989 0.0723684 0.0004689947
#> 
#> $loadings
#>               PC1          PC2          PC3
#> mpg  -0.038118199  0.009184847  0.982070847
#> cyl   0.012035150 -0.003372487 -0.063483942
#> disp  0.899568146  0.435372320  0.031442656
#> hp    0.434784387 -0.899307303  0.025093049
#> drat -0.002660077 -0.003900205  0.039724928
#> wt    0.006239405  0.004861023 -0.084910258
#> qsec -0.006671270  0.025011743 -0.071670457
#> vs   -0.002729474  0.002198425  0.004203328
#> am   -0.001962644 -0.005793760  0.054806391
#> gear -0.002604768 -0.011272462  0.048524372
#> carb  0.005766010 -0.027779208 -0.102897231
#> 
#> $est_dim
#> [1] 3
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
