
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cojam

[![Build
Status](https://travis-ci.com/simisc/cojam.svg?branch=master)](https://travis-ci.com/simisc/cojam)

Multi-variant colocalisation with summary genetic association data.
Merge
[`R2BGLiMS`](https://github.com/pjnewcombe/R2BGLiMS "R2BGLiMS")`::JAM()`
results for two distinct traits, for posterior inference about shared
genetic variants associated with both traits.

## Installation

You can install cojam from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simisc/cojam")
```

## Background

GWAS are generally analyzed one SNP at a time, generating summary data
that obscure the presence of causal variants because many SNPs in
linkage disequilibrium (LD) with the causal SNP will also have
significant associations. JAM (joint analysis of marginal associations,
[`R2BGLiMS` package](https://github.com/pjnewcombe/R2BGLiMS "R2BGLiMS"),
Newcombe et al. 2016) is a scalable fine-mapping algorithm to identify
candidate causal SNPs that best explain the joint pattern of single-SNP
summary associations. The summary associations are conditioned on each
other by estimating their correlation from a reference genotype panel
(e.g. UK Biobank).

Because summary data for different traits are usually available from
distinct studies, `cojam` first models the causal genetic variants for
each trait, then combines these results to draw inferences about the
existence of joint causal variants.

## Combination of JAM model spaces

To do.

## Example

To do.
