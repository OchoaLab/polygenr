
# polygenr

<!-- badges: start -->
<!-- badges: end -->

The goal of `polygenr` is to experiment with raking variants in sparse/penalized regression approaches, particularly `glmnet`.

## Installation

<!-- You can install the released version of polygenr from [CRAN](https://CRAN.R-project.org) with: -->
<!-- install.packages("polygenr") -->

Install from this GitHub repository using:
``` r
library(devtools)
install_github('OchoaLab/polygenr')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(polygenr)

# expected data:
# X: genotype matrix, with loci along rows, individuals along columns
# y: trait vector, length equal to number of individals
# pcs: PC/eigenvector matrix, individuals along rows, dimensions along columns

# fit sparse model while controling for population structure using principal components
obj <- glmnet_pca(X, y, pcs)

# most important component is sparse coefficient matrix:
beta <- obj$beta

# one way of ranking variants, by order of apearance (as penalty factor lambda is decreased):
# returns a vector with length equal to the rows of `beta`
scores <- scores_glmnet( beta )

# another way to rank
# by refitting and calculating ANOVA type-II p-values for each submodel (each column of `beta`)
# returns a matrix of the same dimensions as `beta`
pvals <- anova_glmnet( beta, X, y, pcs )
```

