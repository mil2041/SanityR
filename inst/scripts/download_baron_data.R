# Helper Functions
write_gz_table <- function(x, path) {
    gz <- gzfile(path, open = "w")
    tryCatch(
        write.table(x, gz, col.names = TRUE, row.names = FALSE, quote = FALSE),
        finally = close(gz)
    )
}

# Download UMI counts
# md5:ff15753e2af3cf5271e902f373bc6f0a
# licence: CC-4.0
baron_url <- url("https://zenodo.org/records/4009187/files/Baron_UMI_counts.txt.gz")
baron_counts <- read.table(gzcon(baron_url, text = TRUE), header = TRUE, row.names = 1L)
baron_counts <- as.matrix(baron_counts)

# Cell data
# md5:789dea1233294bf81a0b542564e01af2
# licence: GPL-3 (https://github.com/jmbreda/Sanity)
baron_ct_url <- url("https://raw.githubusercontent.com/jmbreda/Sanity/refs/heads/master/reproducibility/data/Baron_Celltype.txt")
cell_data <- data.frame(
    cell_id  = colnames(baron_counts),
    celltype = scan(baron_ct_url, character(), quiet = TRUE),
    libsize  = colSums(baron_counts)
)
write_gz_table(cell_data, "inst/extdata/baron_cell_data.txt.gz")

# Gene data
gene_data <- data.frame(
    symbol = rownames(baron_counts),
    sum    = rowSums(baron_counts)
)
write_gz_table(gene_data, "inst/extdata/baron_gene_data.txt.gz")
