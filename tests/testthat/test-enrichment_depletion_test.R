context("test-enrichment_depletion_test")

# Read distribution data
distr <- readRDS(system.file("states/distr_data.rds",
  package = "MutationalPatterns"
))
# Set tissue
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))

## Perform the enrichment/depletion test by tissue type.
output <- enrichment_depletion_test(distr, by = tissue)

## Or without specifying the 'by' parameter.
output_pooled <- enrichment_depletion_test(distr)

## Use different cutoffs for p and fdr
output_strictcutoff <- enrichment_depletion_test(distr,
  by = tissue,
  p_cutoffs = 0.000001, fdr_cutoffs = 0.000005
)

# Use multiple cutoffs for p and fdr
output_multistars <- enrichment_depletion_test(distr,
  by = tissue,
  p_cutoffs = c(0.05, 0.01, 0.00000005),
  fdr_cutoffs = c(0.1, 0.05, 0.00000001)
)
test_that("Output has correct class", {
  expect_true(inherits(output, c("data.frame")))
  expect_true(inherits(output_pooled, c("data.frame")))
  expect_true(inherits(output_strictcutoff, c("data.frame")))
  expect_true(inherits(output_multistars, c("data.frame")))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(15, 13))
  expect_equal(dim(output_pooled), c(5, 13))
  expect_equal(dim(output_strictcutoff), c(15, 13))
  expect_equal(dim(output_multistars), c(15, 13))
})

test_that("Number significant is correct", {
  expect_equal(sum(output$significant == "*"), 15)
  expect_equal(sum(output$significant_fdr == "*"), 15)
  expect_equal(sum(output_pooled$significant == "*"), 5)
  expect_equal(sum(output_pooled$significant_fdr == "*"), 5)
  expect_equal(sum(output_strictcutoff$significant == "*"), 9)
  expect_equal(sum(output_strictcutoff$significant_fdr == "*"), 9)
  expect_equal(sum(output_multistars$significant == "***"), 8)
  expect_equal(sum(output_multistars$significant_fdr == "**"), 9)
})
