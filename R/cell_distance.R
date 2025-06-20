#' Calculate the Sanity distance between samples
#'
#' Calculates the expected squared Euclidean distance between two cells using a
#' hierarchical model that shrinks noisy gene differences toward zero.
#'
#' @param x A \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment}
#' object which stores the results of the Sanity analysis.
#' @param assay The name of the assay containing the log normalized counts matrix.
#' @param assay.sd The name of the assay containing the standard deviation of the log-normalized counts
#' @param gene_sd The name of the column in the `rowData(x)` that contains the standard deviation of the gene log-fold change.
#' @param gene_mu The name of the column in the `rowData(x)` that contains the mean log activity of the genes.
#' @param mu_sd The name of the column in the `rowData(x)` that contains the standard deviation of the mean log activity of the genes.
#' @param snr_cutoff A numeric value indicating the minimum signal-to-noise ratio (SNR) to consider a gene.
#' @param nbin Number of bins to use when calculating prior variance of the true distance.
#' @param subset.row A vector of row indices or logical vector indicating which rows to use.
#' @param BPPARAM A BiocParallelParam object specifying the parallelization strategy.
#'
#' @return A \linkS4class{dist} object containing the expected pairwise distances between cells.
#'
#' @details
#'
#' ## Distance Calculation
#'
#' The method calculates the expected squared Euclidean distance between two cells,
#' adjusting for uncertainty in gene expression estimates. For each gene \eqn{g},
#' the contribution to the squared distance between cells \eqn{c} and \eqn{c'} is:
#'
#' \deqn{\langle \Delta_g^2 \rangle = x_g^2 f_g^2(\alpha) + \eta_g^2 f_g(\alpha)}
#'
#' where:
#' - \eqn{x_g = \delta_{gc} - \delta_{gc'}} (observed difference in Sanity's estimates)
#' - \eqn{\eta_g^2 = \epsilon_{gc}^2 + \epsilon_{gc'}^2} (combined error variance)
#' - \eqn{f_g(\alpha) = \alpha v_g/(\alpha v_g + \eta_g^2)} (shrinkage factor)
#'
#' The shrinkage factor balances the observed gene expression differences
#' \eqn{x_g} against their measurement uncertainty \eqn{\eta_g}. For genes with
#' high-confidence estimates (\eqn{\eta_g \rightarrow 0}), it preserves the
#' observed differences while for noisy genes (\eqn{\eta_g \gg 0}), it shrinks
#' the result towards the common expected biological variation inferred from the
#' data (\eqn{\alpha v_g}).
#'
#' The function returns the square root of the expected squared distance
#' \deqn{\langle d \rangle = \sqrt{\sum_g \langle \Delta_g^2 \rangle}}
#'
#' ## Hyperparameter \eqn{\alpha}
#'
#' The key hyperparameter \eqn{\alpha} controls the prior distribution of \eqn{\Delta_g}:
#'
#' \deqn{\Delta_g \sim N(0, \alpha v_g)}
#'
#' Thus:
#'   - \eqn{\alpha = 0}: the 2 cells have identical expression states.
#'   - \eqn{\alpha = 2}: the 2 cells have independent expression states.
#'
#' The function implements numerical integration over \eqn{\alpha} using a grid
#' of `nbin` values to compute the expected value of the squared distance across
#' all possible \eqn{\alpha}.
#'
#' ## Single to Noise Ratio (SNR)
#'
#' *Signal-to-Noise Ratio* (SNR) is defined as the ratio of the variance of
#' log-normalized counts across cells versus the mean variance (i.e. error bars)
#' for each genes.
#'
#' @examples
#' sce <- simulate_branched_random_walk(N_gene = 500, N_path = 10, length_path = 10)
#' sce <- Sanity(sce)  # necessary step before computing distances
#' d <- calculateSanityDistance(sce)
#'
#' # Downstream analysis and visualization
#' hc <- hclust(d, method = "ward.D2")
#' plot(hc)
#'
#' @export
#' @importFrom scuttle .subset2index
#' @importFrom SummarizedExperiment assay rowData
calculateSanityDistance <- function(x,
                                    assay = "logcounts",
                                    assay.sd = "logcounts_sd",
                                    gene_sd = "sanity_activity_sd",
                                    gene_mu = "sanity_log_activity_mean",
                                    mu_sd = "sanity_log_activity_mean_sd",
                                    snr_cutoff = 1, nbin = 400L,
                                    subset.row = NULL, BPPARAM = bpparam()) {

    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x, byrow = TRUE)
        x <- x[subset.row, ]
    }

    delta <- assay(x, assay) - rowData(x)[[gene_mu]]
    epsilon <- assay(x, assay.sd)^2 - rowData(x)[[mu_sd]]^2
    gene_var <- rowData(x)[[gene_sd]]^2

    dmat <- .calculate_sanity_distance(delta, epsilon, gene_var,
                                       snr_cutoff, nbin, BPPARAM)
    attr(dmat, "call") <- match.call()
    return(dmat)
}

#' @importFrom MatrixGenerics rowVars rowMeans2 colMeans2
#' @importFrom BiocParallel bpparam bpmapply
#' @importFrom utils combn
.calculate_sanity_distance <- function(delta, epsilon, gene_var,
                                       snr_cutoff = 1, nbin = 401L,
                                       BPPARAM = bpparam()) {
    # Check validity of inputs
    stopifnot("Dimension Mismatch between mean and errorbars of delta" =
                  all(dim(delta) == dim(epsilon)))
    stopifnot("Delta errorbars cannot be negative" = all(epsilon >= 0))
    stopifnot("Number of rows (genes) in delta matrix must match length of gene_var" =
                  nrow(delta) == length(gene_var))
    stopifnot("Gene variance cannot be negative" = all(gene_var >= 0))
    stopifnot("SNR cutoff cannot be negative" = snr_cutoff >= 0)

    alpha <- seq(0, 2, length.out = nbin + 1L)[-1L]  # skip 0

    if (snr_cutoff > 0) {
        snr <- rowVars(delta) / rowMeans2(epsilon)
        keep <- snr >= snr_cutoff
        delta <- delta[keep, ]
        epsilon <- epsilon[keep, ]
        gene_var <- gene_var[keep]
    }

    ## Rescale (eq. 62-63)
    factor <- gene_var / pmax(gene_var - epsilon, 1e-6)
    delta <- delta * factor
    epsilon <- epsilon * factor
    # clean up before major computational step
    rm(factor)
    gc(verbose = FALSE)

    Ncells <- ncol(delta)
    cell_pairs <- combn(Ncells, 2L)

    get_sanity_distance <- function(i, j) {
        x <- (delta[, i] - delta[, j])^2    # observed distance squared
        eta <- epsilon[, i] + epsilon[, j]  # variance of observed distance assuming independent δ_{gc}

        # Compute variances of true distance Δ_g
        prior_var <- outer(gene_var, alpha, "*")  # a*v_g (GxK)
        posterior_var <- prior_var + eta

        # Compute posterior likelihood  marginalized over alpha (eq. 67 is Supp)
        lik <- -.5 * colMeans2(x / posterior_var + log(posterior_var))  # loglik of 0-centered Gaussian
        lik <- exp(lik - max(lik))  # for numerical stability
        lik <- lik / sum(lik)  # assuming uniform prior over alpha

        # Compute mean distance square from eq 72
        shrinkage <- prior_var / posterior_var  # f_g(a): eq 71
        d2 <- shrinkage * (shrinkage * x + eta)  # (f*x)^2 + f*eta^2
        as.numeric(colMeans2(d2) %*% lik)  # Expected distance squared
    }

    d2 <- bpmapply(get_sanity_distance, cell_pairs[1L, ], cell_pairs[2L, ],
                   SIMPLIFY = TRUE, USE.NAMES = FALSE, BPPARAM = BPPARAM)

    d <- structure(sqrt(d2),
                   Size = Ncells, Labels = colnames(delta),
                   Diag = FALSE, Upper = FALSE,
                   method = "sanity", call = match.call())
    class(d) <- "dist"
    return(d)
}
