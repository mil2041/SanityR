#' Simulate SingleCellExperiment Datasets with Independent or Branched Gene Expression Patterns
#'
#' These functions generate synthetic single-cell RNA-seq datasets based on the
#' Bayesian Model defined in Sanity.
#'
#' @param cell_size Optional vector of real or simulated total UMI counts per cell.
#'   If `NULL`, defaults to values from the *Baron et al.* study.
#' @param gene_size Optional vector of real or simulated total UMI counts per gene.
#'   If `NULL`, defaults to values from the *Baron et al.* study.
#' @param N_cell Integer. Number of cells to simulate.
#' (For `simulate_branched_random_walk` is equal to `N_path * length_path`).
#' If `NULL` inferred from `cell_size`.
#' @param N_gene Integer. Number of genes to simulate. If `NULL`, inferred from `gene_size`.
#' @param seed Integer. Random seed for reproducibility.
#' @param ltq_var_rate Rate parameter for the exponential distribution used to
#' simulate per-gene variance (default: `0.5`).
#' @param N_path (Only for `simulate_branched_random_walk`) Number of branching paths (default: `149`).
#' @param length_path (Only for `simulate_branched_random_walk`) Number of steps (cells) per path (default: `13`).
#'
#' @return A `SingleCellExperiment` object containing:
#' \itemize{
#'   \item \code{assays$counts}: Simulated UMI count matrix.
#'   \item \code{assays$logFC}: Simulated log fold-changes for each gene-cell pair.
#'   \item \code{rowData}: Gene-level metadata including `ltq_mean` and `ltq_var`.
#'   \item \code{colData}: Cell-level metadata including `predecesor` for `simulated_branched_random_walk`.
#' }
#'
#' @details
#'
#' - `simulate_independent_cells`: gene expression values are generated
#' independently for each cell. This results in uncorrelated expression patterns
#' across the dataset.
#'
#' - `simulate_branched_random_walk`: gene expression follows a **branched random walk**
#' through gene expression space, producing correlated gene expression patterns
#' between cells that reflect pseudo-temporal differentiation trajectories.
#'
#' @references
#' A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure
#' Baron, Maayan et al. *Cell Systems*, Volume 3, Issue 4, 346 - 360.e4 \url{https://doi.org/10.1016/j.cels.2016.08.011}
#'
#' @examples
#' # Simulate dataset with independent gene expression
#' sce_indep <- simulate_independent_cells(N_cell = 100, N_gene = 50)
#'
#' # Simulate dataset with a branched random walk trajectory
#' sce_branch <- simulate_branched_random_walk(N_path = 20, length_path = 5, N_gene = 50)
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom SummarizedExperiment "colData<-"
#' @importFrom S4Vectors DataFrame
#' @name simulate_sce
NULL

# Internal helper for argument resolution
# Loads default values if not provided, infers number of cells/genes, and samples sizes
.resolve_simulation_inputs <- function(cell_size, gene_size, N_cell, N_gene, length_path = NULL) {
    if (is.null(cell_size)) {
        path <- system.file("extdata", "baron_cell_data.txt.gz", package = "SanityR")
        cell_size <- read.table(gzfile(path), header = TRUE, stringsAsFactors = FALSE)$libsize
    }

    if (is.null(gene_size)) {
        path <- system.file("extdata", "baron_gene_data.txt.gz", package = "SanityR")
        gene_size <- read.table(gzfile(path), header = TRUE, stringsAsFactors = FALSE)$sum
    }

    if (!is.null(N_cell) && !is.null(length_path)) {
        N_cell <- N_cell * length_path
    } else if (is.null(N_cell)) {
        N_cell <- length(cell_size)
    }

    if (is.null(N_gene)) {
        N_gene <- length(gene_size)
    }

    cell_size <- sample(cell_size, N_cell, replace = TRUE)
    gene_size <- sample(gene_size, N_gene, replace = TRUE)

    list(cell_size = cell_size,
         gene_size = gene_size,
         N_cell = N_cell,
         N_gene = N_gene)
}

# Core simulation logic
# Shared function that computes normalized logFC, transcription rates, counts, and builds SCE
.simulate_sce_core <- function(delta, cell_size, gene_size, ltq_var_rate) {
    N_gene <- nrow(delta)
    N_cell <- ncol(delta)

    # Step 1: Simulate LTQ means and variances
    ltq_mean <- log(gene_size / sum(gene_size))  # mu_g
    ltq_var <- rexp(N_gene, rate = ltq_var_rate) # var_g

    # Step 2: Normalize logFC matrix (delta) to match gene-level variance exactly
    scaling_factor <- sqrt(ltq_var / apply(delta, 1L, var))
    delta <- (delta - rowMeans(delta)) * scaling_factor

    # Step 3: Compute transcription rates and simulate counts
    ltq <- ltq_mean + delta - 0.5 * ltq_var  # logNormal mean = exp(mu + var/2)
    tx_rates <- exp(ltq) %*% diag(cell_size)
        counts_matrix <- matrix(rpois(N_gene * N_cell, tx_rates), nrow = N_gene)

    rownames(counts_matrix) <- rownames(delta) <- paste0("Gene_", seq_len(N_gene))
    colnames(counts_matrix) <- colnames(delta) <- paste0("Cell_", seq_len(N_cell))

    row_data <- DataFrame(ltq_mean = ltq_mean, ltq_var = ltq_var,
                          row.names = rownames(counts_matrix))
    col_data <- DataFrame(cell_size = cell_size,
                          row.names = colnames(counts_matrix))

    sce <- SingleCellExperiment(
        assays = list(counts = counts_matrix, logFC = delta),
        rowData = row_data,
        colData = col_data
    )

    keep <- rowSums(counts(sce)) > 0
    sce <- sce[keep, ]
    return(sce)
}

#' @export
simulate_independent_cells <- function(cell_size = NULL, gene_size = NULL,
                                       N_cell = NULL, N_gene = NULL,
                                       ltq_var_rate = 0.5) {
    args <- .resolve_simulation_inputs(cell_size, gene_size, N_cell, N_gene)
    with(args, {
        # Generate random uncorrelated gene expression matrix
        delta <- matrix(rnorm(N_gene * N_cell), nrow = N_gene)
        .simulate_sce_core(delta, cell_size, gene_size, ltq_var_rate)
    })
}

#' @export
simulate_branched_random_walk <- function(cell_size = NULL, gene_size = NULL,
                                          N_gene = NULL, ltq_var_rate = 0.5,
                                          N_path = 149L, length_path = 13L) {
    args <- .resolve_simulation_inputs(cell_size, gene_size, N_path, N_gene, length_path)
    with(args, {
        # Generate correlated expression through a branched random walk
        delta <- matrix(0, nrow = N_gene, ncol = N_cell)
        predecesor <- integer(N_cell)
        cell <- 1L

        for (k in seq_len(N_path)) {
            # Pick parent node
            if (k == 1L) { # Root
                predecesor[cell] <- 0L
                d_0 <- rep(0, N_gene)
            } else {
                predecesor[cell] <- sample(cell - 1L, 1L)
                d_0 <- delta[, predecesor[cell]]
            }

            # Random walk
            for (i in seq_len(length_path)) {
                if (i == 1L) {
                    delta[, cell] <- d_0 + rnorm(N_gene)
                } else {
                    predecesor[cell] <- cell - 1L
                    delta[, cell] <- delta[, predecesor[cell]] + rnorm(N_gene)
                }
                cell <- cell + 1L
            }
        }

        sce <- .simulate_sce_core(delta, cell_size, gene_size, ltq_var_rate)
        colData(sce)[["predecesor"]] <- predecesor

        return(sce)
    })
}
