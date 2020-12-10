
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cojam

[![R build
status](https://github.com/simisc/cojam/workflows/R-CMD-check/badge.svg)](https://github.com/simisc/cojam/actions)
[![DOI](https://zenodo.org/badge/191818616.svg)](https://zenodo.org/badge/latestdoi/191818616)
[![Licence](https://img.shields.io/github/license/simisc/cojam)](https://github.com/simisc/cojam/blob/master/LICENSE)
[![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg)](https://www.tidyverse.org/lifecycle/#dormant)

Multi-variant colocalisation with summary genetic association data.
Merge [JAM](https://github.com/pjnewcombe/R2BGLiMS "R2BGLiMS package")
results for two distinct traits, for posterior inference about shared
genetic variants associated with both traits.

## Installation

You can install cojam from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simisc/cojam")
```

## Background

GWAS are generally analysed one SNP at a time, generating summary data
that obscure the presence of causal variants because many SNPs in
linkage disequilibrium (LD) with the causal SNP will also have
significant associations.
[JAM](https://github.com/pjnewcombe/R2BGLiMS "R2BGLiMS package") is a
scalable fine-mapping algorithm to identify candidate causal SNPs that
best explain the joint pattern of single-SNP summary associations
([Newcombe et al. 2016 *Genetic
Epidemiology*](https://doi.org/10.1002/gepi.21953 "JAM paper 1"),
[Newcombe et al. 2018 *Genetic
Epidemiology*](https://doi.org/10.1002/gepi.22245 "JAM paper 2")). The
summary associations are conditioned on each other by estimating their
correlation from a reference genotype panel (e.g. UK Biobank).

Because summary data for different traits are usually available from
distinct studies, cojam first models the causal genetic variants for
each trait, then combines these results to draw inferences about the
existence of joint causal variants. Building on the colocalisation
framework used in [coloc](https://github.com/chr1swallace/coloc "coloc")
([Giambartolomei et al. 2014 *PLOS
Genetics*](https://doi.org/10.1371/journal.pgen.1004383 "coloc paper 1"),
[Wallace 2020 *PLOS
Genetics*](https://doi.org/10.1371/journal.pgen.1008720 "coloc paper 2")),
cojam assigns joint models for the two traits to one of five mutually
exclusive hypotheses. The two hypotheses of interest are:

  - H<sub>3</sub>: associations with both traits, distinct SNPs
    (colocalisation)
  - H<sub>4</sub>: associations with both traits, including at least one
    shared SNP

In cojam, support for H<sub>4</sub> over H<sub>3</sub> is quantified by
the Bayes Factor (likelihood ratio), *BF* = p(D|H<sub>4</sub>) /
p(D|H<sub>3</sub>) = PosteriorOdds\[H<sub>4</sub>:H<sub>3</sub>\] /
PriorOdds\[H<sub>4</sub>:H<sub>3</sub>\]. The posterior odds are
estimated from JAM’s independent rjMCMC posterior samples; the prior
odds are calculated combinatorially, also assuming independence between
traits, and taking into account the (potentially different) priors on
the proportion of causal variants in each JAM model. The effects of the
independence assumption on the prior and posterior odds cancel out in
the *BF*.

When two traits both have genetic drivers in the same region, this
suggests some relationship between the traits. For a chosen prior odds
of H<sub>4</sub> against H<sub>3</sub> (reflecting the dependence
between colocalised traits), multiplying by *BF* provides the resulting
posterior odds.

**Caution**: The current version of cojam does not include methods for
assessing whether merged chains have efficiently searched the joint
model space, or for visualising multi-SNP results. Both JAM models
should be thoroughly checked before passing them to cojam, using methods
provided in
[R2BGLiMS](https://github.com/pjnewcombe/R2BGLiMS "R2BGLiMS package").

## Example

To do.
