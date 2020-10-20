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

tags <- function(genotypes_1, genotypes_2 = NULL, threshold = 0.9) {

    r2a <- r2b <- genotypes_1 %>%
        rbind(genotypes_2) %>%
        cov()

    hc <- hclust(as.dist(1 - r2), "single")
    clusters <- cutree(hc, h = 1 - threshold ^ 2)
    r2a[outer(clusters, clusters, "!=")] <- NA
    r2b[outer(clusters, clusters, "==")] <- NA

    enframe(clusters, name = "snp", value = "group") %>%
        mutate(r2_in = apply(r2a, 1, mean, na.rm = TRUE),
               r2_out = apply(r2b, 1, mean, na.rm = TRUE)) %>%
        group_by(group) %>%
        mutate(group_size = n(),
               rank = order(order(-r2_in, r2_out))) %>%
        filter(snp %in% vars1,
               snp %in% vars2) %>%
        filter(rank == min(rank)) %>%
        "$"(snp)
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
        genotype_matrix <- na.omit(genotype_matrix)
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
