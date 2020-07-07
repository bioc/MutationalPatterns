context("test-plot_contribution_heatmap")


# Read in nmf results
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
  package = "MutationalPatterns"
))
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")

# Plot with clustering.
output <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples = TRUE, cluster_sigs = TRUE)

# Define signature and sample order for plotting.
sig_order <- c("Signature B", "Signature A")
sample_order <- c(
  "colon1", "colon2", "colon3", "intestine1", "intestine2",
  "intestine3", "liver3", "liver2", "liver1"
)
output_supplied_order <- plot_contribution_heatmap(nmf_res$contribution,
  cluster_samples = FALSE,
  sig_order = sig_order, sample_order = sample_order
)

## Contribution heatmap with text values
output_text <- plot_contribution_heatmap(nmf_res$contribution, plot_values = TRUE)

# Read in signature refitting results
snv_refit <- readRDS(system.file("states/strict_snv_refit.rds",
  package = "MutationalPatterns"
))
output_refit <- plot_contribution_heatmap(snv_refit$contribution, cluster_samples = TRUE, cluster_sigs = TRUE)


test_that("Output has correct class", {
  expect_true(inherits(output, "gg"))
  expect_true(inherits(output_supplied_order, "gg"))
  expect_true(inherits(output_text, "gg"))
  expect_true(inherits(output_refit, "gg"))
})
