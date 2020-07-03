context("test-plot_cosine_heatmap")

# Get mut_mat
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Get signatures
signatures <- get_known_signatures()


# Calculate the cosine similarity between each signature and each 96 mutational profile
cos_matrix <- cos_sim_matrix(mut_mat, signatures)

# Plot the cosine similarity between each signature and each sample with hierarchical
# clustering of samples and signatures.
output <- plot_cosine_heatmap(cos_matrix, cluster_rows = TRUE, cluster_cols = TRUE, method = "complete")

# In the above example, clustering is performed on the similarities of the samples with
# the signatures. It's also possible to cluster the signatures and samples on their (96) profile.
hclust_cosmic <- cluster_signatures(signatures, method = "average")
cosmic_order <- colnames(signatures)[hclust_cosmic$order]
hclust_samples <- cluster_signatures(mut_mat, method = "average")
sample_order <- colnames(mut_mat)[hclust_samples$order]
# Plot the cosine heatmap using this given signature order.
output_supplied_order <- plot_cosine_heatmap(cos_matrix,
  cluster_rows = FALSE, cluster_cols = FALSE, row_order = sample_order,
  col_order = cosmic_order, method = "complete"
)

# You can also plot the similarity of samples with eachother
cos_matrix <- cos_sim_matrix(mut_mat, mut_mat)
output_inner <- plot_cosine_heatmap(cos_matrix, cluster_rows = TRUE, cluster_cols = TRUE, method = "complete")

# You can also include test
output_text <- plot_cosine_heatmap(cos_matrix, cluster_rows = TRUE, cluster_cols = TRUE, plot_values = TRUE)

test_that("Output has correct class", {
  expect_true(inherits(output, "gg"))
  expect_true(inherits(output_supplied_order, "gg"))
  expect_true(inherits(output_inner, "gg"))
  expect_true(inherits(output_text, "gg"))
})
