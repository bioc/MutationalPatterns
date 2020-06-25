context("test-plot_signature_strand_bias")

# Load strand data
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

# Load nmf results
nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
  package = "MutationalPatterns"
))

## Provide column names for the plot.
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")

output <- plot_signature_strand_bias(nmf_res_strand$signatures)

# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
})
