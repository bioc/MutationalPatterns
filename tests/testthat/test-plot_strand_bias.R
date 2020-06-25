context("test-plot_strand_bias")

# Read stranded mut_mat
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

## Perform the strand bias test.
strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
strand_bias <- strand_bias_test(strand_counts)

## Plot the strand bias.
output <- plot_strand_bias(strand_bias)

# Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
  package = "MutationalPatterns"
))
strand_counts <- strand_occurrences(mut_mat_repli, by = tissue)
strand_bias <- strand_bias_test(strand_counts)
output_repli <- plot_strand_bias(strand_bias)

## Test with p instead of fdr
output_pval <- plot_strand_bias(strand_bias, sig_type = "p")


## Use multiple (max 3) significance cutoffs.
strand_bias_multistars <- strand_bias_test(strand_counts,
  p_cutoffs = c(0.05, 0.01, 0.005),
  fdr_cutoffs = c(0.1, 0.05, 0.01)
)
output_multistars <- plot_strand_bias(strand_bias_multistars)


test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_repli, c("gg")))
  expect_true(inherits(output_pval, c("gg")))
  expect_true(inherits(output_multistars, c("gg")))
})
