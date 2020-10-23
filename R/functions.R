#' Wrapper around \code{\link[R2BGLiMS]{JAM}}
#'
#' Constructs a collated list of data that can be passed to \code{\link[R2BGLiMS]{JAM}}
#'   via \code{do.call}, or as an argument to \code{\link{cojam}}. The prior
#'   proportion of causal covariates is treated as unknown and given a Beta(a, b)
#'   hyper-prior.
#'
#' @param marginal_beta Vector of marginal effect estimates to re-analyse with
#'   JAM under multivariate models. For GWAS summaries, these or \emph{log-ORs}.
#' @param snp_names SNP identifiers (e.g. RSID), in the same order as \code{marginal_beta}.
#' @param ref_genotypes Reference genotype matrix used by JAM to impute the SNP-SNP
#'   correlations. Individual's genotype must be coded as a numeric risk allele
#'   count 0/1/2. Non-integer values reflecting imputation may be given. NB: The
#'   risk allele coding MUST correspond to that used in marginal.betas. Must be
#'   positive definite, with SNP identifiers in the column names.
#' @param n The size of the dataset in which the summary statistics
#'   \code{marginal.betas} were calculated
#' @param trait_variance Estimate of the trait (outcome) variance.
#' @param binary_outcome Is the trait (outcome) binary? Should be \code{TRUE} for
#'   GWAS or \code{FALSE} for QTL studies.
#' @param marginal_beta_se Only required if the trait (outcome) is binary:
#'   standard errors of the log-ORs.
#' @param bb_prior_a Parameter of Beta(a, b) prior on proportion of causal
#'   covariates, default 1.
#' @param bb_prior_b Parameter of Beta(a, b) prior on proportion of causal
#'   covariates. Default is the number of SNPs that are common to both \code{snp_names}
#'   and \code{colnames(ref_genotypes)}. Higher values of \code{bb_prior_b}
#'   relative to \code{bb_prior_a} will encourage greater sparsity. Use
#'   \code{\link[R2BGLiMS]{GetBetaBinomParams}} for suggested beta-binomial parameters.
#' @param ... Other arguments to \code{\link[R2BGLiMS]{JAM}}.
#' @return A collated list of data that can be passed to \code{\link[R2BGLiMS]{JAM}}
#'   via \code{do.call}.
#' @export
jam_wrap <- function(marginal_beta,
                     snp_names,
                     ref_genotypes,
                     n,
                     trait_variance,
                     binary_outcome = FALSE,
                     marginal_beta_se = NULL,
                     bb_prior_a = 1,
                     bb_prior_b = NULL,
                     ...) {


    variables <- intersect(snp_names, colnames(ref_genotypes)) # Add message if this results in loss of variables...
    if (is.null(bb_prior_b)) {
        bb_prior_b <- length(variables)
    }

    names(marginal_beta) <- snp_names
    marginal_beta <- marginal_beta[variables]
    ref_genotypes <- ref_genotypes[, variables]

    extra_arguments <- NULL # for continuous outcome, default inverse Gamma prior, a = b = 0.01

    if (binary_outcome) {

        is.null(marginal_beta_se) &&
            stop("For binary outcomes, supply marginal_beta_se: std errors of the log-ORs")

        names(marginal_beta_se) <- snp_names
        marginal_beta_se <- marginal_beta_se[variables]

        # NB: extra arguments must be named!
        extra_arguments <- list(GaussianResidualVarianceInvGammaPrior_a = 2,
                                GaussianResidualVarianceInvGammaPrior_b = trait_variance)

        z_score <- (marginal_beta / marginal_beta_se) / sqrt(n)
        marginal_beta <- z_score * sqrt(trait_variance) / apply(ref_genotypes, 2, stats::sd)
    }

    list(marginal.betas = marginal_beta,
         X.ref = ref_genotypes,
         n = n,
         trait.variance = trait_variance,
         model.space.prior = list(
             a = bb_prior_a,
             b = bb_prior_b,
             Variables = variables
         ),
         extra.arguments = extra_arguments,
         ...)
}

jam_models <- function(jam_res) {
    jam_res@mcmc.output %>%
        tibble::as_tibble() %>%
        dplyr::select(-alpha, -LogLikelihood) %>%
        dplyr::group_by_all() %>%
        dplyr::summarise(posterior = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(posterior)) %>%
        dplyr::mutate(model_size = rowSums(.[, jam_res@model.space.prior$Variables]),
                      model_rank = dplyr::row_number())
}

#' Beta-Binomial probabilities for specific models
#'
#' This is not the prior on a particular dimension (that would required the
#'   binomial co-efficient). For use in calculating the implied prior (under
#'   independence between traits) in \code{cojam}.
#'
#' @param jam_res A \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object,
#'   from running \code{\link[R2BGLiMS]{JAM}}.
#' @return Vector of the same length P (number of SNPs in \code{jam_res}),
#'   containing prior probabilities for models of dimensions 0 to (P-1).
model_size_priors <- function(jam_res) {
    a <- jam_res$model.space.prior$a
    b <- jam_res$model.space.prior$b
    n <- length(jam_res$model.space.prior$Variables)
    k <- 0:(n - 1) # dimension of specific model
    beta(k + a, n - k + b) / beta(a, b)
}

implied_prior <- function(jam_res_1, jam_res_2) {

    n_snps <- length(jam_res_1$model.space.prior$Variables)

    model_priors1 <- model_size_priors(jam_res_1)
    model_priors2 <- model_size_priors(jam_res_2)

    # p_H0 <- model_priors_1[1] * model_priors_2[1]
    # p_H1 <- model_priors_2[1] - p_H0
    # p_H2 <- model_priors_1[1] - p_H0
    p_H0_U_H1_U_H2 <- model_priors_1[1] + model_priors_2[1] - model_priors_1[1] * model_priors_2[1]

    p_H3 <- expand.grid(i = 1:(n_snps - 1),
                        j = 1:(n_snps - 1)) %>%
        dplyr::filter(i + j <= n_snps) %>%
        dplyr::mutate(prior_1 = model_priors_1[i + 1],
                      prior_2 = model_priors_2[j + 1],
                      h3_combos = choose(n_snps, i) * choose(n_snps - i, j),
                      h3_prior = prior_1 * prior_2 * h3_combos) %>%
        dplyr::pull(h3_prior) %>%
        sum()

    p_H4 <- 1 - p_H0_U_H1_U_H2 - p_H3

    c(H012 = p_H0_U_H1_U_H2, H3 = p_H3, H4 = p_H3)
}

#' Merge two \code{\link[R2BGLiMS]{JAM}} models
#'
#' Combines samples from the posterior model space of two \code{\link[R2BGLiMS]{JAM}}
#' models, allowing posterior inference about colocalisation, i.e. about the
#' existence of one or more shared variants that are associated with both outcome
#' traits.
#'
#' @param jam_res_1 A \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object,
#'   from running \code{\link[R2BGLiMS]{JAM}} for trait 1.
#' @param jam_res_2 A \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object,
#'   from running \code{\link[R2BGLiMS]{JAM}} for trait 2.
#' @param prior_odds_43 Optional vector of prior odds of a shared variant ("H4") versus
#'   distinct variants ("H3"). Can also be added later using \code{\link{cojam_odds}}.
#' @return A \code{list} containing the following elements:
#'   \describe{
#'       \item{\code{bayes_factor}:}{Bayes Factor for H4 versus H3.}
#'       \item{\code{posterior_overlap}:}{Posterior probability of H4 U H3.}
#'       \item{\code{posterior_odds_43}:}{If \code{prior_odds_43} are provided,
#'       the posterior odds of H4 versus H3.}
#'       \item{\code{models_grid}:}{All samples configurations (joint
#'         \code{\link[R2BGLiMS]{JAM}} models) with posterior counts.}
#'     }
#' @export
cojam <- function(jam_res_1, jam_res_2, prior_odds_43 = NA) {

    vars <- jam_res_1$model.space.prior$Variables

    identical(vars, jam_res_2$model.space.prior$Variables) ||
        stop("JAM calls must include identical variables. Use tag() to find a suitable set of SNPs.")
    (jam_res_1@enumerate.up.to.dim == 0 & jam_res_2@enumerate.up.to.dim == 0) ||
        stop("Not yet implemented for enumerated JAM models.")
    (jam_res_1@n.covariate.blocks.for.jam == 1 & jam_res_2@n.covariate.blocks.for.jam == 1) ||
        stop("Not yet implemented for multiple covariate blocks.")

    jam_tab1 <- jam_models(jam_res_1)
    jam_tab2 <- jam_models(jam_res_2)

    models_grid <- tcrossprod(as.matrix(jam_tab1[, vars]),
                              as.matrix(jam_tab2[, vars])) %>%
        reshape2::melt(varnames = c("model_rank_1", "model_rank_2"),
                       value.name = "num_coloc") %>%
        tibble::as_tibble() %>%
        dplyr::full_join(jam_tab1, by = "model_rank_1") %>%
        dplyr::full_join(jam_tab2, by = "model_rank_2", suffix = c("_1", "_2")) %>%
        dplyr::mutate(
            hypoth = dplyr::case_when(
                # model_size_1 == 0 & model_size_2 == 0 ~ "H0",
                # model_size_2 == 0 ~ "H1",
                # model_size_1 == 0 ~ "H2",
                model_size_1 == 0 | model_size_2 == 0 ~ "H012",
                num_coloc == 0 ~ "H3",
                TRUE ~ "H4"),
            posterior = posterior_1 * posterior_2)

    posterior_counts <- tibble::tibble(hypoth = c("H012", "H3", "H4"), by = "hypoth") %>%
        dplyr::full_join(models_grid) %>%
        dplyr::group_by(hypoth) %>%
        dplyr::summarise(posterior = sum(posterior, na.rm = TRUE), .groups = "drop") %>%
        dplyr::pull(posterior, name = hypoth)

    prior_probs_implied <- implied_prior(jam_res_1, jam_res_2)
    prior_odds_34_implied <- prior_probs_implied["H3"] / prior_probs_implied["H4"]

    posterior_overlap <- sum(posterior_counts[c("H3", "H4")]) / sum(posterior_counts)
    posterior_odds_43_implied <- posterior_counts["H4"] / posterior_counts["H3"]

    bayes_factor <- posterior_odds_43_implied * prior_odds_34_implied

    posterior_odds_43 <- tibble::tibble(prior = prior_odds_43) %>%
        dplyr::mutate(posterior = prior * bayes_factor)

    list(bayes_factor = bayes_factor,
         posterior_overlap = posterior_overlap,
         posterior_odds_43 = posterior_odds_43,
         models_grid = models_grid)
}

#' Augment \code{\link{cojam}} results with other prior choices
#'
#' Calculate posterior odds of a shared variant ("H4") versus distinct variants ("H3")
#'
#' @param cojam_res Results from \code{\link{cojam}}
#' @param prior_odds_43 Vector of prior odds of H4 versus H3.
#' @return Same as \code{\link{cojam}}, a \code{list} containing the following elements:
#'   \describe{
#'       \item{\code{bayes_factor}:}{Bayes Factor for H4 versus H3.}
#'       \item{\code{posterior_overlap}:}{Posterior probability of H4 U H3.}
#'       \item{\code{posterior_odds_43}:}{Posterior odds of H4 versus H3, including
#'       those already in \code{cojam_res} (if any) and new rows for new \code{prior_odds_43}.}
#'     }
#' @export
cojam_odds <- function(cojam_res, prior_odds_43) {
    add_posterior_odds_43 <- tibble::tibble(prior = prior_odds_43) %>%
        dplyr::mutate(posterior = prior * cojam_res$bayes_factor)
    cojam_res$posterior_odds_43 <- cojam_res$posterior_odds_43 %>%
        dplyr::bind_rows(add_posterior_odds_43) %>%
        dplyr::distinct()
    cojam_res
}

plot_shared_variants <- function(cojam_res) {
    # Manhattan plots for trait 1 and trait 2
    # Manhattan plot for PP of shared variant(s) *given* H4
    # Have partial code for this somewhere...
    # ggplot2::theme_classic() + ggplot2::theme(panel.border = element_rect(fill = NA))
}

#' Pipe
#'
#' Imported from \code{\link[magrittr]{pipe}}
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

#' Plot a correlation matrix
#'
#' Pretty correlation matrix with `ggplot`
#'
#' @param M Correlation matrix
#' @export
plot_cormat <- function(M) {

    isSymmetric(M) || stop("Correlation matrix must be symmetric")

    M %>%
        reshape2::melt() %>%
        ggplot2::ggplot() +
        ggplot2::geom_raster(ggplot2::aes(Var1, stats::reorder(Var2, dplyr::desc(Var2)), fill = value)) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 90),
                       panel.border = element_rect(fill = NA)) +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
        ggplot2::coord_fixed()
}

#' Convert DNA string to its complement
#'
#' Convert string of nucleotide codes to its complement
#'
#' @param x Character string of nucleotide code
#' @return Character string of the complement of \code{x}
#' @export
complement_DNA <- function(x) {

    letters_x <- unlist(stringr::str_split(x, ""))
    letters_DNA <- "ACGTMRWSYKVHDBN-" # from Biostrings::DNA_ALPHABET

    all(stringr::str_detect(letters_DNA, letters_x)) ||
        stop(sprintf("Acceptable input characters are %s", letters_DNA))

    chartr(letters_DNA, "TGCAKYWSRMBDHVN-", x)
}

#' Hierarchical pruning of genotype matrix
#'
#' Find SNPs that "tag" groups of correlated SNPs by cutting a hierarchical
#'   clustering tree at the specified threshold.
#'
#' @param genotypes_1 Genotype matrix, with SNP identifiers as column names.
#' @param genotypes_2 Optional second genotype matrix, useful to jointly select
#'   tag SNPs for use in \code{cojam}.
#' @param threshold Correlation (height) at which to cut the hierarchical
#'   clustering tree.
#' @return Character vector of tag SNP identifiers.
#' @export
prune_tags <- function(genotypes_1, genotypes_2 = NULL, threshold = 0.9) {

    if (!is.null(genotypes_2)) {
        vars <- intersect(colnames(genotypes_1), colnames(genotypes_2))
        genotypes_1 <- genotypes_1[, vars] %>%
            rbind(genotypes_2[, vars])
    }

    r2a <- r2b <- stats::cov(genotypes_1)

    hc <- stats::hclust(stats::as.dist(1 - r2a), "single")
    clusters <- stats::cutree(hc, h = 1 - threshold ^ 2)
    r2a[outer(clusters, clusters, "!=")] <- NA
    r2b[outer(clusters, clusters, "==")] <- NA

    tibble::enframe(clusters, name = "snp", value = "group") %>%
        dplyr::mutate(r2_within = apply(r2a, 1, mean, na.rm = TRUE),
                      r2_between = apply(r2b, 1, mean, na.rm = TRUE)) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(group_size = dplyr::n(),
                      rank = order(order(-r2_within, r2_between))) %>%
        dplyr::filter(rank == min(rank)) %>%
        dplyr::pull(snp)
}

#### Legacy functions ####

logsum <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

prune_qr <- function(genotypes) {

    if (any(is.na(genotypes))) {
        n <- nrow(genotypes)
        genotypes <- stats::na.omit(genotypes)
        message(sprintf(
            "Removing %i rows containing missing values (prune_qr).",
            n - nrow(genotypes)
        ))
    }

    gm_qr <- qr(genotypes)
    genotypes[, gm_qr$pivot[seq_len(gm_qr$rank)]]
}

prune_pairwise <- function(genotypes, threshold = 0.9) {

    genotypes <- prune_qr(genotypes)

    cormat <- abs(stats::cor(genotypes))
    cormat[cormat > abs(threshold)] <- NA
    diag(cormat) <- 1

    while(any(is.na(cormat))) {

        remove_snp <- tibble::tibble(
            label = rownames(cormat),
            num_NAs = apply(is.na(cormat), 1, sum),
            avg_cor = apply(cormat, 1, mean, na.rm = TRUE)
        ) %>%
            dplyr::filter(num_NAs == max(num_NAs)) %>%
            dplyr::filter(avg_cor == max(avg_cor)) %>%
            dplyr::pull(label)

        length(remove_snp) == 1 || stop("Tied SNPs in in prune_genotypes?")

        cormat <- cormat[rownames(cormat) != remove_snp,
                         rownames(cormat) != remove_snp,
                         drop = FALSE]
    }

    genotypes[, rownames(cormat), drop = FALSE]
}

# to appease R CMD check
utils::globalVariables(c("LogLikelihood", "Var1", "Var2", "alpha", "avg_cor", "group",
                         "h3_combos", "h3_prior", "hypoth", "i", "j", "label",
                         "model_priors_1", "model_priors_2", "num_NAs", "posterior",
                         "posterior_1", "posterior_2", "prior", "prior_1", "prior_2",
                         "r2_between", "r2_within", "snp", "value", "."))
