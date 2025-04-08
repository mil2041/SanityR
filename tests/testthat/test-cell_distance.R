library(testthat)
library(SummarizedExperiment)
library(BiocParallel)

set.seed(2025)
sce <- simulate_branched_random_walk(N_gene = 12, N_path = 3, length_path = 5)
sce <- Sanity(sce)

test_that("calculateSanityDistance returns a valid dist object", {
    dist_obj <- calculateSanityDistance(
        sce,
        assay = "logcounts",
        assay.sd = "logcounts_sd",
        gene_sd = "sanity_activity_sd",
        gene_mu = "sanity_log_activity_mean",
        mu_sd = "sanity_log_activity_mean_sd",
        snr_cutoff = 0,
        nbin = 10L,
        BPPARAM = SerialParam()
    )

    expect_true(inherits(dist_obj, "dist"))
    expect_equal(attr(dist_obj, "Size"), ncol(sce))
    expect_equal(length(dist_obj), choose(ncol(sce), 2))
})

test_that("calculateSanityDistance works with subset.row", {
    # Select only the first 6 genes.
    subset_rows <- 1:6
    dist_obj_subset <- calculateSanityDistance(
        sce,
        assay = "logcounts",
        assay.sd = "logcounts_sd",
        gene_sd = "sanity_activity_sd",
        gene_mu = "sanity_log_activity_mean",
        mu_sd = "sanity_log_activity_mean_sd",
        snr_cutoff = 0.0,
        nbin = 10L,
        subset.row = subset_rows,
        BPPARAM = SerialParam()
    )

    expect_true(inherits(dist_obj_subset, "dist"))
    expect_equal(attr(dist_obj_subset, "Size"), ncol(sce))
    expect_equal(length(dist_obj_subset), choose(ncol(sce), 2))

    # Compare with full data calculation to ensure using subset.row alters the result.
    dist_obj_all <- calculateSanityDistance(
        sce,
        assay = "logcounts",
        assay.sd = "logcounts_sd",
        gene_sd = "sanity_activity_sd",
        gene_mu = "sanity_log_activity_mean",
        mu_sd = "sanity_log_activity_mean_sd",
        snr_cutoff = 0.0,
        nbin = 10L,
        BPPARAM = bpparam()
    )
    expect_false(identical(as.numeric(dist_obj_subset), as.numeric(dist_obj_all)))
})
