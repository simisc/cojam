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
        marginal_beta <- z_score * sqrt(trait_variance) / apply(ref_genotypes, 2, sd)
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

jam_models <- function(jamres) {

    msp <- jamres@model.space.priors[[1]]

    if (jamres@enumerate.up.to.dim > 0) {
        stop("Not yet implemented for enumerated JAM models.")
        # names(msp$Variables) <- msp$Variables
        # model_probs <- jamres@enumerated.posterior.inference$model.probs
        # restab <- map_dfr(msp$Variables, function(v) {
        #     stringr::str_detect(names(model_probs), v)
        # })
        # restab$posterior <- model_probs
    }

    jamres@mcmc.output %>%
        tibble::as_tibble() %>%
        dplyr::select(-alpha, logLike = LogLikelihood) %>%
        dplyr::group_by_all() %>%
        dplyr::summarise(posterior = n(), .groups = "drop") %>%
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
    a <- jam_arg$model.space.prior$a
    b <- jam_arg$model.space.prior$b
    n <- length(jam_arg$model.space.prior$Variables)
    k <- 0:(n - 1)
    beta(k + a, n - k + b) / beta(a, b)
}

implied_prior <- function(jam_arg1, jam_arg2) {

    if (!identical(jam_arg1$model.space.prior$Variables, jam_arg2$model.space.prior$Variables)) {
        stop("JAM calls must include identical variables. Use tag() to find a suitable set of SNPs.")
    }

    n_snps <- length(jam_arg1$model.space.prior$Variables)

    model_priors1 <- model_size_priors(jam_arg1)
    model_priors2 <- model_size_priors(jam_arg2)

    p_H0 <- model_priors_1[1] * model_priors_2[1]
    p_H1 <- model_priors_2[1] - p_H0
    p_H2 <- model_priors_1[1] - p_H0
    p_H3_U_H4 <- 1 - p_H2 - p_H1 - p_H0

    p_H3 <- expand.grid(i = 1:(n_snps - 1),
                        j = 1:(n_snps - 1)) %>%
        filter(i + j <= n_snps) %>%
        mutate(prior_1 = model_priors_1[i + 1],
               prior_2 = model_priors_2[j + 1],
               h3_combos = choose(n_snps, i) * choose(n_snps - i, j),
               h3_prior = prior_1 * prior_2 * h3_combos) %>%
        pull(h3_prior) %>%
        sum()

    p_H4 <- p_H3_U_H4 - p_H3

    c(H0 = p_H0, H1 = p_H1, H2 = p_H2, H3 = p_H3, H4 = p_H3)
}

# prior_odds_43 can be a vector of prior odds: p(H4) / p(H3)
cojam <- function(jam_arg1, jam_arg2, prior_odds_43 = NA) {

    if (!identical(jam_arg1$model.space.prior$Variables, jam_arg2$model.space.prior$Variables)) {
        stop("JAM calls must include identical variables. Use tag() to find a suitable set of SNPs.")
    }

    jam_res1 <- do.call(JAM, jam_arg1)
    jam_res2 <- do.call(JAM, jam_arg2)

    jam_tab1 <- jam_models(jam_res1)
    jam_tab2 <- jam_models(jam_res2)

    posterior_counts <- tcrossprod(as.matrix(jam_tab1[, v1]),
                                   as.matrix(jam_tab2[, v1])) %>%
        reshape2::melt(varnames = c("model_rank_1", "model_rank_2"),
                       value.name = "num_coloc") %>%
        tibble::as_tibble() %>%
        dplyr::full_join(jam_tab1, by = "model_rank_1") %>%
        dplyr::full_join(jam_tab2, by = "model_rank_2", suffix = c("_1", "_2")) %>%
        dplyr::mutate(
            hypoth = dplyr::case_when(
                model_size_1 == 0 & model_size_2 == 0 ~ "H0",
                model_size_2 == 0 ~ "H1",
                model_size_1 == 0 ~ "H2",
                num_coloc == 0 ~ "H3",
                TRUE ~ "H4"),
            posterior_indep = posterior_1 * posterior_2
        ) %>%
        full_join(
            tibble::tibble(hypoth = paste0("H", 0:4), by = "hypoth")
        ) %>%
        group_by(hypoth) %>%
        summarise(posterior_indep = sum(posterior_indep, na.rm = TRUE), .groups = "drop") %>%
        pull(posterior_indep, name = hypoth)

    prior_probs_implied <- implied_prior(jam_arg1, jam_arg2)
    prior_odds_34_implied <- prior_probs["H3"] / prior_probs["H4"]

    posterior_overlap <- sum(posterior_counts[c("H3", "H4")]) / sum(posterior_counts)
    posterior_odds_43_implied <- posterior_counts["H4"] / posterior_counts["H3"]

    bayes_factor <- posterior_odds_43_implied * prior_odds_34_implied

    posterior_odds_43 <- tibble::tibble(prior = prior_odds_43) %>%
        dplyr::mutate(posterior = prior * bayes_factor)

    list(bayes_factor = bayes_factor,
         posterior_overlap = posterior_overlap,
         posterior_odds_43 = posterior_odds_43)
}
