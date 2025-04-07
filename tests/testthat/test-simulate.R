library(testthat)
library(SingleCellExperiment)

set.seed(42)

test_that("simulate_independent_cells returns a valid SCE", {
    sce <- simulate_independent_cells(N_cell = 20, N_gene = 15)

    # Check that the returned object is a SingleCellExperiment
    expect_s4_class(sce, "SingleCellExperiment")

    # Check that required assays exist and have correct dimensions
    expect_true(all(c("counts", "logFC") %in% assayNames(sce)))
    expect_true(all(rowSums(counts(sce)) > 0))
    expect_lte(nrow(sce), 15)
    expect_equal(ncol(sce), 20)

    # Check that rowData and colData contain the expected columns
    expect_true(all(c("ltq_mean", "ltq_var") %in% colnames(rowData(sce))))
    expect_true(all("cell_size" %in% colnames(colData(sce))))
})

test_that("simulate_branched_random_walk returns a valid SCE with lineage graph", {
    sce <- simulate_branched_random_walk(N_gene = 10, N_path = 3, length_path = 5)

    # Check that the returned object is a SingleCellExperiment
    expect_s4_class(sce, "SingleCellExperiment")

    # Check that required assays exist
    expect_true(all(c("counts", "logFC") %in% assayNames(sce)))
    expect_true(all(rowSums(counts(sce)) > 0))
    expect_lte(nrow(sce), 15)
    expect_equal(ncol(sce), 3 * 5)

    # Check that colData contains neighbor info
    expect_true(all(c("cell_size", "predecesor") %in% colnames(colData(sce))))
})
