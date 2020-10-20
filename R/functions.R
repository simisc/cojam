#' Pipe
#'
#' Imported from \code{\link[magrittr]{pipe}}
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

prune_qr <- function(genotype_matrix) {
    if (any(is.na(genotype_matrix))) {
        n <- nrow(genotype_matrix)
        genotype_matrix <- na.omit(genotype_matrix)
        message(sprintf(
            "Removing %i rows containing missing values (prune_qr).",
            n - nrow(genotype_matrix)
        ))
    }

    gm_qr <- qr(genotype_matrix)
    genotype_matrix[, gm_qr$pivot[seq_len(gm_qr$rank)]]
}

beta_binom_specific_model <- function(k, n, a, b) {
    # Calculates a Beta-Binomial prior probability for a SPECIFIC model. From Bottolo et al.
    # This is not the prior on a particular dimension (that would required the binomial co-efficient).
    # k dimension of specific model
    # n total number of covariates
    # a beta-binomial hyper parameter a
    # b beta-binomial hyper parameter b
    # returns Probability
    beta(k + a, n - k + b) / beta(a, b)
}

logsum <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

jam_model_table <- function(jamres) {

    msp <- jamres@model.space.priors[[1]]

    if (jamres@enumerate.up.to.dim == 0) {

        restab <- jamres@mcmc.output %>%
            tibble::as_tibble() %>%
            dplyr::select(-alpha, logLike = LogLikelihood) %>%
            dplyr::group_by_all() %>%
            dplyr::summarise(postProb = n()) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(postProb = postProb / sum(postProb))

    } else {
        stop("Not yet implemented for enumerated JAM models (cojam priors are wrong)")
        names(msp$Variables) <- msp$Variables
        model_probs <- jamres@enumerated.posterior.inference$model.probs
        restab <- map_dfr(msp$Variables, function(v) {
            stringr::str_detect(names(model_probs), v)
            })
        restab$postProb <- model_probs
    }

    restab %>%
        dplyr::arrange(dplyr::desc(postProb)) %>%
        dplyr::mutate(modelDim = rowSums(.[, msp$Variables]),
                      modelRank = dplyr::row_number())
}

subset_jam_args <- function(jam_args, vars) {

    if(!all(vars %in% jam_args$model.space.prior$Variables)) {
        stop("Not all vars are in the existing JAM calls (subset_jam_args).")
    }

    jam_args$marginal.betas <- jam_args$marginal.betas[vars]
    jam_args$X.ref <- jam_args$X.ref[, vars]
    jam_args$model.space.prior$b <- jam_args$model.space.prior$b * length(vars) / length(jam_args$model.space.prior$Variables)
    jam_args$model.space.prior$Variables <- vars

    jam_args
}

hypoth_priors <- function(lambda1 = 1, lambda2 = 1, n_snps = 1000) {

    pk1 <- lambda1 / (1 + lambda1)                    # p(k1 = 0), null model 1
    pk2 <- lambda2 / (1 + lambda2)                    # p(k2 = 0), null model 2

    odds_implied_prior <- pk1 * pk2 * n_snps
    # p(H3) / p(H4) under independence
    # Based on linear APPROXIMATION in size P, not good for P<25 ??

    pr3_3u4_indep <- odds_implied_prior / (1 + odds_implied_prior)
    # p(H3 | {H3 U H4}) under independence

    pr012 <- c(pk1 * pk2,                               # p(H0)
               (1 - pk1) * pk2,                         # p(H1)
               pk1 * (1 - pk2))                         # p(H2)
    pr3u4 <- (1 - pk1) * (1 - pk2)                      # p(H3 U H4)
    pr34 <- pr3u4 * c(pr3_3u4_indep, 1 - pr3_3u4_indep) # p(H3), p(H4)

    list(
        hypoth_priors = tibble::tibble(hypoth = 0:4, prior_indep = c(pr012, pr34)),
        colocalisation_priors = c(odds_colocalisation = odds_implied_prior, # p(H3) / p(H4)
                                  probability_overlap = pr3u4)              # p(H3 U H4)
    )
}

cojam <- function(jam_arg1, jam_arg2, prior_odds = NULL) {

    if(jam_arg1$model.space.prior$a != 1 | jam_arg2$model.space.prior$a != 1) {
        stop("Model space priors of both JAM calls must be of the form BetaBin(1, b).")
    }

    v1 <- jam_arg1$model.space.prior$Variables
    v2 <- jam_arg2$model.space.prior$Variables

    if(!identical(v1, v2)) {

        v1 <- intersect(v1, v2)

        warning(sprintf("Both JAM calls must use same SNPs\nKeeping intersection: %i SNPs", length(vars)))

        jam_arg1 <- subset_jam_args(jam_arg1, v1)
        jam_arg2 <- subset_jam_args(jam_arg2, v1)
    }

    lambda1 <- jam_arg1$model.space.prior$b / length(jam_arg1$model.space.prior$Variables)
    lambda2 <- jam_arg2$model.space.prior$b / length(jam_arg2$model.space.prior$Variables)

    jam_res1 <- do.call(JAM, jam_arg1)
    jam_res2 <- do.call(JAM, jam_arg2)

    jam_tab1 <- jam_model_table(jam_res1)
    jam_tab2 <- jam_model_table(jam_res2)

    X <- tcrossprod(as.matrix(jam_tab1[, v1]),
                    as.matrix(jam_tab2[, v1]))

    names(jam_tab1) <- sprintf("%s_1", names(jam_tab1))
    names(jam_tab2) <- sprintf("%s_2", names(jam_tab2))

    res_grid <- reshape2::melt(X,
                               varnames = c("modelRank_1", "modelRank_2"),
                               value.name = "num_coloc") %>%
        tibble::as_tibble() %>%
        dplyr::left_join(jam_tab1, by = "modelRank_1") %>%
        dplyr::left_join(jam_tab2, by = "modelRank_2") %>%
        dplyr::mutate(
            hypoth = dplyr::case_when(
                modelDim_1 == 0 & modelDim_2 == 0 ~ 0,
                modelDim_2 == 0 ~ 1,
                modelDim_1 == 0 ~ 2,
                num_coloc == 0 ~ 3,
                TRUE ~ 4),
            postprob_indep = postProb_1 * postProb_2
        )

    implied_prior_independent <- hypoth_priors(lambda1 = lambda1,
                                               lambda2 = lambda2,
                                               n_snps = length(v1))

    hypotheses_summary <- implied_prior_independent$hypoth_priors %>%
        dplyr::group_by_all() %>%
        dplyr::full_join(res_grid, by = "hypoth") %>%
        dplyr::summarise(postprob_indep = sum(postprob_indep, na.rm = TRUE)) %>%
        dplyr::ungroup()

    # Need implied prior here after all...
    postprob34 <- c(H4 = hypotheses_summary$postprob_indep[hypotheses_summary$hypoth == 4],
                    H3 = hypotheses_summary$postprob_indep[hypotheses_summary$hypoth == 3])
    posterior_odds <- postprob34["H4"] / postprob34["H3"]
    bf <- posterior_odds * implied_prior_independent$colocalisation_priors$odds_colocalisation

    res <- list(bayes_factor = bf,
                posterior_probability_overlap = sum(postprob34),
                posterior_probabilities = NULL,       # tibble with column for each prior
                posterior_odds_colocalisation = NULL, # tibble(prior_odds, posterior_odds)
                jam_results = list(jam1 = jam_res1,
                                   jam2 = jam_res2,
                                   grid = res_grid),
                implied_prior_independent = implied_prior_independent,
                pars = list(lambda1 = lambda1, lambda2 = lambda2,
                            variables = v1, n_snps = length(v1)))

    if (!is.null(prior_odds)) {
        res <- cojam_postprobs(res, prior_odds)
    }

    res
}

cojam_postprobs <- function(cojam_res, prior_odds) {

    odds_colocalisation = tibble(prior = prior_odds) %>%
        mutate(posterior = prior * cojam_res$bayes_factor)


    # Post odds and prior odds should be in the same direction!

    # Rearrange res to hold all independent results in one element, then
    # add a new element for each prior. Can then add elements (using another)
    # call to cojm_postprobs) if not already included?

    # res_summ <- priors_tab %>%
    #     dplyr::full_join(res_grid, by = "hypoth") %>%
    #     dplyr::mutate(
    #         postprob_indep = postProb_1 * postProb_2,
    #         postprob_joint = postprob_indep * weight,
    #
    #         ## Workaround: is this correct? Or should I be calculating posterior ##
    #         ## probabilities directly from likelihoods and (reweighted) priors?  ##
    #         postprob_joint = postprob_joint / sum(postprob_joint, na.rm = TRUE)
    #     ) %>%
    #     dplyr::group_by(hypoth, weight, prior_indep, prior_joint) %>%
    #     dplyr::summarise(
    #         postprob_indep = sum(postprob_indep, na.rm = TRUE),
    #         postprob_joint = sum(postprob_joint, na.rm = TRUE)
    #     ) %>%
    #     ungroup()

}

prune_genotypes <- function(genotypes, threshold = 0.8) {

    warning("Obsolete: use hiearchical version instead (prune_genotypes).")

    genotypes <- prune_qr(genotypes)

    cormat <- abs(cor(genotypes))
    cormat[cormat > abs(threshold)] <- NA
    diag(cormat) <- 1

    while(any(is.na(cormat))) {

        remove_snp <- tibble(
            label = rownames(cormat),
            num_NAs = apply(is.na(cormat), 1, sum),
            avg_cor = apply(cormat, 1, mean, na.rm = TRUE)
            ) %>%
            filter(num_NAs == max(num_NAs)) %>%
            filter(avg_cor == max(avg_cor)) %>%
            "$"(label)

        length(remove_snp) == 1 || stop("Tied SNPs in in prune_genotypes?")

        cormat <- cormat[rownames(cormat) != remove_snp,
                         rownames(cormat) != remove_snp,
                         drop = FALSE]
    }

    genotypes[, rownames(cormat), drop = FALSE]
}

plot_cormat <- function(M) {
    isSymmetric(M) || stop("Correlation matrix must be symmetric")

    M %>%
        reshape2::melt() %>%
        ggplot2::ggplot() +
        # ggplot2::geom_raster(aes(Var1, Var2, fill = value)) +
        ggplot2::geom_raster(aes(Var1, reorder(Var2, desc(Var2)), fill = value)) +
        ggplot2::theme(axis.title = element_blank(),
                       axis.text.x = element_text(angle = 90)) +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
        ggplot2::coord_fixed()
}

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

complement_DNA <- function(x) {
    letters_x <- unlist(stringr::str_split(x, ""))
    letters_DNA <- "ACGTMRWSYKVHDBN-" # from Biostrings::DNA_ALPHABET

    if (!all(stringr::str_detect(letters_DNA, letters_x))) {
        stop(sprintf("Acceptable input characters are %s", letters_DNA))
    }

    chartr(letters_DNA, "TGCAKYWSRMBDHVN-", x)
}
