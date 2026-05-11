
<!-- README.md is generated from README.Rmd. Please edit that file -->

<style>
p {
  text-align: justify;
}
</style>

# nlkernel

<!-- badges: start -->

<!-- badges: end -->

This package provides a unified framework for genomic prediction using
kernel-based methods within the RKHS paradigm implemented in BGLR. It
includes functions for fitting and evaluating Laplacian, Gaussian,
Bessel and Polynomial kernels. It also allows combining multiple kernels
in a single model. Still, these kernels can be used in Principal
Component Analysis (KPCA), with automatic selection of principal
components based on explained variance thresholds. Finally, it allows
combining multiple kernels in a single model.

The package supports flexible hyperparameter tuning to all kernels.

For comparison purposes, it is possible to implement GBLUP and PCA in
the same prediction structure.

The models are available for both single-environment and
multi-environment trials. For the latter, frameworks are available
considering (1) main effects only, and (2) both main and interaction
effects.

The results returned the accuracy, measured as the correlation between
predicted and observed values.

## Installation and execution

Package installation can be done directly by calling:

``` r
pak::pak("wwanessa13/nlkernel")
```

To load the package into memory, use:

``` r
library(nlkernel)
```

To see all functions avaliable, run:

``` r
ls("package:nlkernel")
```

More details about functions are show running:

``` r
?function_name
```

## Single environment

To evaluate model performance, K-fold cross-validation is used.

### Functions names

- gblup
- laplacian
- gaussian
- bessel
- polynomial
- combinations
- pca
- pca_laplacian
- pca_gaussian
- pca_bessel
- pca_polynomial

### Example with settings default

``` r
SNPs <- read.table("gen.txt")

y <- read.table("phenotypes.txt") %>% pull(yield)

results <- gaussian(SNPs, y)
print(results)
```

## Multi environment

To evaluate model performance, different cross-validation (CV) schemes
can be applied:

- CV1: prediction performance of unobserved genotypes in observed
  environments.

- CV2: prediction performance of genotypes observed in only a subset of
  environments.

- CV0: prediction performance of observed genotypes in unobserved
  environments.

### Y = E + G + e

The model treats environments as fixed effects and genotypes as random
effects, following a main-effects genomic framework.

#### Functions name

- env_g_gblup
- env_g_laplacian
- env_g_gaussian
- env_g_bessel
- env_g_polynomial
- env_g_combinations
- env_g_pca
- env_g_pca_laplacian
- env_g_pca_gaussian
- env_g_pca_bessel
- env_g_pca_polynomial

#### Example: Default settings

``` r
SNPs <- read.csv("gen.csv")

phen <- read.csv("Phenotypes.csv")

y <- phen$yield

IDs <- phen$genotype

env <- phen$Env

results <- env_g_gaussian(SNPs, y, IDs, env, CV = "CV1")
print(results)
```

### Y = E + G + GxE + e

Environment is treated as fixed, while G and G×E are random effects.
This main-and-interaction-effects framework applies the kernel to the G
matrix and extends it to the G×E interaction matrix.

#### Functions name

- env_ge_gblup
- env_ge_laplacian
- env_ge_gaussian
- env_ge_bessel
- env_ge_polynomial
- env_ge_combinations
- env_ge_pca
- env_ge_pca_laplacian
- env_ge_pca_gaussian
- env_ge_pca_bessel
- env_ge_pca_polynomial

#### Example: Default settings

``` r
SNPs <- read.csv("gen.csv")

phen <- read.csv("Phenotypes.csv")

y <- phen$yield

IDs <- phen$genotype

env <- phen$Env

results <- env_ge_gaussian(SNPs, y, IDs, env, CV = "CV1")
print(results)
```

## Attention notes

When you run a function with `file_name = NULL` argument, the result is
saved by default using the model name. However, if you run the same
model for more than one trait, the second output will overwrite the
first one. Therefore, when using the same model for multiple traits, we
recommend renaming the `.xlsx` file before running the next trait.

Alternatively, you can set `save_xlsx = FALSE` in the model function
arguments. Then, after inspecting the results, you can save them in your
preferred format and with your desired file name.
