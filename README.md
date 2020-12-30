
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICBioMark

<!-- badges: start -->
<!-- badges: end -->

Welcome to ICBioMark, based on the paper “Data-driven design of targeted
gene panels forestimating immunotherapy biomarkers”, *Bradley and
Cannings* (currently being written
[here](https://github.com/cobrbra/TargetedPanelEstimation_Paper)).
ICBioMark is designed to implement regularised optimisation methods for
the design and use of targeted gene panels in predicting
immunotherapeutic biomarkers such as Tumour Mutation Burden (TMB) and
Tumour Indel Burden (TIB).

## Installation

<!-- You can install the development version of ICBioMark from [CRAN](https://CRAN.R-project.org) with: -->

You can install the development version of this package from this github
repository (using the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package) with:

``` r
devtools::install_github("cobrbra/ICBioMark")
```

## Example

Upon installation we can load the package.

``` r
library(ICBioMark)
## basic example code
```

To demonstrate the typical workflow for using ICBioMark, we play around
with a small and friendly example dataset. This is saved into the
package, but just comes from the data simulation function
`{r eval = FALSE} generate_maf_data()`.

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
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!
