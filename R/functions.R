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
    genotype_matrix[, gm_qr$pivot[seq(gm_qr$rank)]]
}

BetaBinomialProbabilitySpecificModel <- function(k, n, a, b) {
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

    tab <- jamres@mcmc.output %>%
        tibble::as_tibble() %>%
        dplyr::select(-alpha, logLike = LogLikelihood) %>%
        dplyr::group_by_all() %>%
        dplyr::summarise(postProb = n()) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(postProb)) %>%
        dplyr::mutate(
            postProb = postProb / sum(postProb),
            # express as probabilties
            modelDim = rowSums(.[, msp$Variables]),
            logPrior = log(
                BetaBinomialProbabilitySpecificModel(
                    k = modelDim,
                    n = length(msp$Variables),
                    a = msp$a,
                    b = msp$b
                )
            ),
            logPost  = logPrior + logLike,

            # Normalise posterior:
            logPost  = logPost - logsum(logPost),
            estProb  = exp(logPost),
            modelRank = dplyr::row_number()
        )

    attr(tab, "vars") <- msp$Variables
    tab
}

jam_model_grid <- function(jamres1,
                           jamres2,
                           priors = list(lambda1 = 1, lambda2 = 1, odds_colocalisation = 100),
                           vars = NULL) {

    tab1 <- jam_model_table(jamres1)
    tab2 <- jam_model_table(jamres2)

    if (is.null(vars)) {
        v1 <- attr(tab1, "vars")
        v2 <- attr(tab2, "vars")
        vars <- intersect(v1, v2)
    }

    for (i in 1:length(priors)) {
        assign(names(priors)[i], priors[[i]])
    }

    odds_implied_prior <- 1000 * (lambda1 / (1 + lambda1)) * (lambda2 / (1 + lambda2))
    weight3 <- (odds_colocalisation / (odds_colocalisation + 1)) * ((odds_implied_prior + 1) / odds_implied_prior)
    weight4 <- (odds_implied_prior + 1) / (odds_colocalisation + 1)

    cat(
        sprintf(
            "Implied odds: %.2f\nWeight 3: %.2f\nWeight 4: %.2f\n",
            odds_implied_prior, weight3, weight4
        )
    )

    X <- tcrossprod(as.matrix(tab1[, vars]),
                    as.matrix(tab2[, vars]))

    names(tab1) <- sprintf("%s_1", names(tab1))
    names(tab2) <- sprintf("%s_2", names(tab2))

    restab <-
        reshape2::melt(X,
                       varnames = c("modelRank_1", "modelRank_2"),
                       value.name = "num_coloc") %>%
        tibble::as_tibble() %>%
        dplyr::left_join(tab1, by = "modelRank_1") %>%
        dplyr::left_join(tab2, by = "modelRank_2") %>%
        dplyr::mutate(
            hypoth = dplyr::case_when(
                modelDim_1 == 0 & modelDim_2 == 0 ~ 0,
                modelDim_2 == 0 ~ 1,
                modelDim_1 == 0 ~ 2,
                num_coloc == 0 ~ 3,
                TRUE ~ 4
            ),
            reweighting_factor = dplyr::case_when(hypoth == 4 ~ weight4, # 251 / 101,
                                                  hypoth == 3 ~ weight3, # 2510 / 2525,
                                                  TRUE ~ 1),

            logLike_joint  = logLike_1 + logLike_2,

            # Assuming independence (no reweighting):
            logPrior_indep = logPrior_1 + logPrior_2,
            logPost_indep  = logPrior_indep + logLike_joint,
            # Normalise posterior:
            logPost_indep  = logPost_indep - logsum(logPost_indep),
            estProb_indep  = exp(logPost_indep),

            # Reweighting the prior depending on coloc hypotheses:
            logPrior_joint = logPrior_indep + log(reweighting_factor)
        ) %>%
        dplyr::mutate(
            logPost_joint  = logPrior_joint + logLike_joint,
            # Normalise posterior:
            logPost_joint  = logPost_joint - logsum(logPost_joint),
            estProb_joint  = exp(logPost_joint)
        )
    restab

}

cojam <- function(jam_arg1, jam_arg2, prior_odds = 100, vars = NULL) {

    if(jam_arg1$model.space.prior$a != 1 | jam_arg2$model.space.prior$a != 1) {
        stop("Model space priors of both JAM calls must be of the form BetaBin(1, b).")
    }

    jam_res1 <- do.call(JAM, jam_arg1)
    jam_res2 <- do.call(JAM, jam_arg2)

    lambda1 <- jam_arg1$model.space.prior$b / length(jam_arg1$model.space.prior$Variables)
    lambda2 <- jam_arg2$model.space.prior$b / length(jam_arg2$model.space.prior$Variables)

    res <- jam_model_grid(jam_res1, jam_res2, vars = vars,
                          priors = list(lambda1 = lambda1, lambda2 = lambda2, odds_colocalisation = prior_odds))

    #################################################################
    ### Table is wrong: only takes prior sum over visited models! ###
    ### Should calculate the prior, not take it from the MCMC...  ###
    #################################################################

    tibble::tibble(hypoth = 0:4) %>%
        dplyr::full_join(res) %>%
        dplyr::group_by(hypoth) %>%
        dplyr::summarise(
            # prior_indep = sum(exp(logPrior_indep), na.rm = TRUE),
            # prior_joint = sum(exp(logPrior_joint), na.rm = TRUE),
            postprob_indep = sum(estProb_indep, na.rm = TRUE),
            postprob_joint = sum(estProb_joint, na.rm = TRUE)
        )
}

prune_genotypes <- function(genotypes, threshold = 0.975) {

   genotypes <- prune_qr(genotypes) # First remove colinear SNPs

    cormat <- abs(cor(genotypes, use = "pairwise.complete.obs"))
    cormat[cormat > abs(threshold)] <- NA
    diag(cormat) <- 1
    remove_idx <- logical(nrow(cormat))

    while(any(is.na(cormat))) {

        remove_snp <- tibble(
            num_NAs = apply(is.na(cormat), 1, sum, na.rm = TRUE),
            avg_cor = apply(cormat, 1, mean, na.rm = TRUE)
        ) %>%
            mutate(index = row_number()) %>%
            filter(num_NAs == max(num_NAs)) %>%
            filter(avg_cor == max(avg_cor)) %>%
            "$"(index)

        cormat <- cormat[-remove_snp, -remove_snp]
        remove_idx[remove_snp] <- TRUE
    }

    genotypes[, !remove_idx, drop = FALSE]
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

make_extra_args_InvGammaPrior <- function(n_cases,
                                          n_total,
                                          n_prior = 1,
                                          version_2 = TRUE) {

    # Called this list extra_args_InvGammaPrior in clean-for-jam.R
    # Pass this list to JAM via the argument 'extra.arguments='
    # New version Paul sent by email...

    p_cases <- n_cases / n_total
    var_cases <- p_cases * (1 - p_cases)

    if (version_2) {
        res <- list(
            GaussianResidualVarianceInvGammaPrior_a = 2,
            GaussianResidualVarianceInvGammaPrior_b = var_cases
        )
    } else {
        b <- var_cases * (1 + (n_prior - 1) / 2)
        res <- list(
            GaussianResidualVarianceInvGammaPrior_a = b / var_cases + 1,
            GaussianResidualVarianceInvGammaPrior_b = b
        )
    }
    return(res)
}

jam_wrap <- function(marginal.betas, # transformed/linear
                     X.ref,
                     snp_names,
                     n,
                     prior_lambda = 1,
                     trait_variance = NULL,
                     inv_gamma_prior = NULL,
                     thinning_interval = NULL,
                     prior_b = prior_lambda * length(snp_names),
                     list_args = FALSE) {

    if (is.null(names(marginal.betas))) {
        message("Assuming that marginal.betas and snp_names are in the same order.")
        names(marginal.betas) <- snp_names
    } else if (all(snp_names %in% names(marginal.betas))) {
        marginal.betas <- marginal.betas[snp_names]
    } else {
        stop("Not all snp_names are in names(marginal.betas)")
    }

    if (all(snp_names %in% colnames(X.ref))) {
        X.ref = na.omit(X.ref[, snp_names])
    } else {
        stop("Not all snp_names are in names(X.ref)")
    }

    if (prior_b != prior_lambda * length(snp_names)) {
        message("Default prior_b was overridden: prior_lambda will be ignored.")
    }

    jam_args <- list(
        marginal.betas = marginal.betas,
        X.ref = X.ref,
        n = n,
        trait.variance = trait_variance,
        model.space.prior = list(
            a = 1,
            b = prior_b,
            Variables = snp_names
        ),
        extra.arguments = inv_gamma_prior,
        thinning.interval = thinning_interval
    )

    if (list_args) {
        return(jam_args)
    } else {
        return(do.call(JAM, jam_args))
    }
}

jam_plot <- function(jam_model,
                     snp_names,
                     univariate_pvalues,
                     true_signal = NULL,
                     colours = FALSE) {
    # true_signal can be:
    # NULL      SNP with lowest p-value is chosen
    # character SNP(s) named are chosen
    # logical   Directly specifies choice of SNP

    data <-
        tibble::tibble(SNP = snp_names, pvalue = univariate_pvalues) %>%
        dplyr::mutate(neglogp = -log10(pvalue),
                      index = row_number())

    if (is.null(true_signal)) {
        data$signal = data$neglogp == max(data$neglogp, na.rm = TRUE)
    } else if (is.character(true_signal)) {
        data$signal = data$SNP %in% true_signal
    } else if (is.logical(true_signal) &
               length(true_signal) == nrow(data)) {
        data$signal = true_signal
    } else {
        stop("true_signal must be NULL, character or logical.")
    }

    post_prob <- jam_model@posterior.summary.table[, "PostProb"] %>%
        na.omit() %>% # drops alpha, LogLikelihood, ModelSizePartition1 to only keep SNPs
        tibble::tibble(post_prob = ., SNP = names(.)) %>%
        dplyr::inner_join(data) %>%
        dplyr::mutate(
            jam_snp = post_prob == max(post_prob),
            col_code = as.numeric(jam_snp) + 2 * as.numeric(signal),
            col_code = recode(
                col_code,
                `0` = "None",
                `1` = "JAM",
                `2` = "Univariate",
                `3` = "Both"
            )
        ) %>%
        tidyr::gather(key, value, post_prob, neglogp)

    if (colours) {
        outplot <-
            ggplot2::ggplot(post_prob, aes(x = index, y = value, col = col_code)) +
            ggplot2::geom_point() +
            ggplot2::facet_wrap(~ key,
                                scales = "free_y",
                                labeller = ggplot2::labeller(
                                    key = c(neglogp = "Univariate -log10(p)",
                                            post_prob = "Marginal posterior probability")
                                ))
    } else {
        outplot <- ggplot2::ggplot(post_prob, aes(x = index, y = value)) +
            ggplot2::geom_point() +
            ggplot2::facet_wrap(~ key,
                                scales = "free_y",
                                labeller = ggplot2::labeller(
                                    key = c(neglogp = "Univariate -log10(p)",
                                            post_prob = "Marginal posterior probability")
                                )) +
            ggplot2::scale_color_discrete(name = "Top SNP")
    }

    outplot
}

complement_DNA <- function(x) {
    letters_x <- unlist(str_split(x, ""))
    letters_DNA <- "ACGTMRWSYKVHDBN-" # Biostrings::DNA_ALPHABET

    if (!all(str_detect(letters_DNA, letters_x))) {
        stop(sprintf("Acceptable input characters are %s", letters_DNA))
    }

    chartr(letters_DNA, "TGCAKYWSRMBDHVN-", x)
}
