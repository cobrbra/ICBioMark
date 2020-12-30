
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

## Example Workflow

Upon installation we can load the package.

``` r
library(ICBioMark)
```

To demonstrate the typical workflow for using ICBioMark, we play around
with a small and friendly example dataset. This is pre-loaded with the
package, but just comes from the data simulation function
`generate_maf_data()`, so you can play around with datasets of different
sizes and shapes.

### Input Data

Our example dataset, called `example_maf_data`, is a list with two
elements: `maf` and `gene_lengths`. These are the two pieces of data
that you’ll always need to use this package, and they look as follows:

-   `maf` is a data frame in the MAF (mutation annotated format) style.
    For a set of sequenced tumour/normal pairs, this means a table with
    a row for every mutation identified, with columns corresponding to
    properties such as the sample ID for the tumour of origin, the gene,
    chromosome and nucelotide location of the mutation, and the type of
    mutation observed. In the real world, MAF datasets often have lots
    of extra information beyond this, but in our small example we’ve
    just included sample, gene and mutation type (it’s all we’ll need!).
    The top five rows look like this:

``` r
 # example_maf_data <- generate_maf_data()
 knitr::kable(head(example_maf_data$maf, 5), row.names = FALSE)
```

| Tumor\_Sample\_Barcode | Hugo\_Symbol | Variant\_Classification |
|:-----------------------|:-------------|:------------------------|
| SAMPLE\_53             | GENE\_5      | Missense\_Mutation      |
| SAMPLE\_27             | GENE\_7      | Missense\_Mutation      |
| SAMPLE\_73             | GENE\_11     | Silent                  |
| SAMPLE\_47             | GENE\_18     | 3’Flank                 |
| SAMPLE\_36             | GENE\_19     | Missense\_Mutation      |

-   `gene_lengths`, another data frame, this time containing the names
    of genes that you’ll want to include in your modelling and their
    length. Gene length is a complex and subtle thing to define - we
    advise using coding length as defined in the
    [Ensembl](https://www.ensembl.org/index.html) database. For this
    example, however, gene lengths are again randomly chosen:

``` r
  knitr::kable(head(example_maf_data$gene_lengths, 5), row.names = FALSE)
```

| Hugo\_Symbol | max\_cds |
|:-------------|---------:|
| GENE\_1      |      961 |
| GENE\_2      |     1009 |
| GENE\_3      |     1011 |
| GENE\_4      |      976 |
| GENE\_5      |     1016 |

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
