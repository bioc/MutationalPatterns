context("test-plot_enrichment_depletion")

# Read distribution data
distr <- readRDS(system.file("states/distr_data.rds",
  package = "MutationalPatterns"
))
# Set tissue
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

## Perform the enrichment/depletion test.
distr_test <- enrichment_depletion_test(distr, by = tissue)
distr_test2 <- enrichment_depletion_test(distr)

## Plot the enrichment/depletion
output <- plot_enrichment_depletion(distr_test)
output_persample <- plot_enrichment_depletion(distr_test2)

## Test with p instead of fdr
output_pval <- plot_enrichment_depletion(distr_test, sig_type = "p")

## Use multiple (max 3) significance cutoffs.
distr_multistars <- enrichment_depletion_test(distr,
  by = tissue,
  p_cutoffs = c(0.05, 0.01, 0.005),
  fdr_cutoffs = c(0.1, 0.05, 0.01)
)
output_multistars <- plot_enrichment_depletion(distr_multistars)


# Perform tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_persample, c("gg")))
  expect_true(inherits(output_pval, c("gg")))
  expect_true(inherits(output_multistars, c("gg")))
})
