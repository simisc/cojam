jam_wrap <- function(marginal_beta, # _original_ betas (log ORs if binary outcome)
                     snp_names,
                     ref_genotypes,
                     n,
                     trait_variance,
                     binary_outcome = FALSE,
                     marginal_beta_se = NULL, # only required if binary outcome
                     prior_lambda = 1) {      # temporary: need to rethink this!

    names(marginal_beta) <- snp_names
    variables <- intersect(snp_names, colnames(ref_genotypes)) # Add message if this results in loss of variables...

    marginal_beta <- marginal_beta[variables]
    ref_genotypes <- ref_genotypes[, variables]

    extra_arguments <- NULL # for continuous outcome, default inverse Gamma prior, a = b = 0.01

    if (binary_outcome) {

        if (is.null(marginal_beta_se)) {
            stop("For binary outcomes, supply marginal_beta_se: std errors of the log-ORs")
        }

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
             a = 1,
             b = prior_lambda * length(variables),
             Variables = variables
         ),
         extra.arguments = extra_arguments)
}

jam_models <- function(jam_res) {

    msp <- jam_res@model.space.priors[[1]]

    if (jam_res@enumerate.up.to.dim > 0) {
        stop("Not yet implemented for enumerated JAM models.")
        # names(msp$Variables) <- msp$Variables
        # model_probs <- jam_res@enumerated.posterior.inference$model.probs
        # restab <- map_dfr(msp$Variables, function(v) {
        #     stringr::str_detect(names(model_probs), v)
        # })
        # restab$posterior <- model_probs
    }

    jam_res@mcmc.output %>%
        tibble::as_tibble() %>%
        dplyr::select(-alpha, -LogLikelihood) %>%
        dplyr::group_by_all() %>%
        dplyr::summarise(posterior = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(posterior)) %>%
        dplyr::mutate(model_size = rowSums(.[, msp$Variables]),
                      model_rank = dplyr::row_number())
}

# Calculates a Beta-Binomial prior probability for a SPECIFIC model. From Bottolo et al.
# This is not the prior on a particular dimension (that would required the binomial co-efficient).
# k dimension of specific model
# n total number of covariates
# a beta-binomial hyper parameter a
# b beta-binomial hyper parameter b
# returns Probability
model_size_priors <- function(jam_args) {
    a <- jam_args$model.space.prior$a
    b <- jam_args$model.space.prior$b
    n <- length(jam_args$model.space.prior$Variables)
    k <- 0:(n - 1)
    beta(k + a, n - k + b) / beta(a, b)
}

implied_prior <- function(jam_args_1, jam_args_2) {

    if (!identical(jam_args_1$model.space.prior$Variables, jam_args_2$model.space.prior$Variables)) {
        stop("JAM calls must include identical variables. Use tag() to find a suitable set of SNPs.")
    }

    n_snps <- length(jam_args_1$model.space.prior$Variables)

    model_priors1 <- model_size_priors(jam_args_1)
    model_priors2 <- model_size_priors(jam_args_2)

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

# prior_odds_43 can be a vector of prior odds: p(H4) / p(H3)
cojam <- function(jam_args_1, jam_args_2, prior_odds_43 = NA) {

    vars <- jam_args_1$model.space.prior$Variables

    if (!identical(vars, jam_args_2$model.space.prior$Variables)) {
        stop("JAM calls must include identical variables. Use tag() to find a suitable set of SNPs.")
    }

    jam_res1 <- do.call(R2BGLiMS::JAM, jam_args_1)
    jam_res2 <- do.call(R2BGLiMS::JAM, jam_args_2)

    jam_tab1 <- jam_models(jam_res1)
    jam_tab2 <- jam_models(jam_res2)

    posterior_counts <- tcrossprod(as.matrix(jam_tab1[, vars]),
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
            posterior = posterior_1 * posterior_2
        ) %>%
        dplyr::full_join(
            tibble::tibble(hypoth = c("H012", "H3", "H4"), by = "hypoth")
        ) %>%
        dplyr::group_by(hypoth) %>%
        dplyr::summarise(posterior = sum(posterior, na.rm = TRUE), .groups = "drop") %>%
        dplyr::pull(posterior, name = hypoth)

    prior_probs_implied <- implied_prior(jam_args_1, jam_args_2)
    prior_odds_34_implied <- prior_probs_implied["H3"] / prior_probs_implied["H4"]

    posterior_overlap <- sum(posterior_counts[c("H3", "H4")]) / sum(posterior_counts)
    posterior_odds_43_implied <- posterior_counts["H4"] / posterior_counts["H3"]

    bayes_factor <- posterior_odds_43_implied * prior_odds_34_implied

    posterior_odds_43 <- tibble::tibble(prior = prior_odds_43) %>%
        dplyr::mutate(posterior = prior * bayes_factor)

    list(bayes_factor = bayes_factor,
         posterior_overlap = posterior_overlap,
         posterior_odds_43 = posterior_odds_43)
}

#### Helpers ####

#' Pipe
#'
#' Imported from \code{\link[magrittr]{pipe}}
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

logsum <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

tags <- function(genotypes_1, genotypes_2, threshold = 0.9) {

    vars <- intersect(colnames(genotypes_1),
                      colnames(genotypes_2))

    r2a <- r2b <- genotypes_1[, vars] %>%
        rbind(genotypes_2[, vars]) %>%
        stats::cov()

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

plot_cormat <- function(M) {
    isSymmetric(M) || stop("Correlation matrix must be symmetric")

    M %>%
        reshape2::melt() %>%
        ggplot2::ggplot() +
        # ggplot2::geom_raster(aes(Var1, Var2, fill = value)) +
        ggplot2::geom_raster(ggplot2::aes(Var1, stats::reorder(Var2, dplyr::desc(Var2)), fill = value)) +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 90)) +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
        ggplot2::coord_fixed()
}

complement_DNA <- function(x) {

    letters_x <- unlist(stringr::str_split(x, ""))
    letters_DNA <- "ACGTMRWSYKVHDBN-" # from Biostrings::DNA_ALPHABET

    if (!all(stringr::str_detect(letters_DNA, letters_x))) {
        stop(sprintf("Acceptable input characters are %s", letters_DNA))
    }

    chartr(letters_DNA, "TGCAKYWSRMBDHVN-", x)
}

prune_qr <- function(genotype_matrix) {

    if (any(is.na(genotype_matrix))) {
        n <- nrow(genotype_matrix)
        genotype_matrix <- stats::na.omit(genotype_matrix)
        message(sprintf(
            "Removing %i rows containing missing values (prune_qr).",
            n - nrow(genotype_matrix)
        ))
    }

    gm_qr <- qr(genotype_matrix)
    genotype_matrix[, gm_qr$pivot[seq_len(gm_qr$rank)]]
}

prune_pairwise <- function(genotypes, threshold = 0.8) {

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
