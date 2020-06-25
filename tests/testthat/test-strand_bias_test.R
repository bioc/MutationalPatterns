context("test-strand_bias_test")

# Load stranded mutation matrix
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

# Set tissue names
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

## Perform the strand bias test.
strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
output <- strand_bias_test(strand_counts)

# Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
  package = "MutationalPatterns"
))
strand_counts_repli <- strand_occurrences(mut_mat_repli, by = tissue)
output_repli <- strand_bias_test(strand_counts_repli)

## Use different cutoffs for p and fdr
output_lenientcutoff <- strand_bias_test(strand_counts, p_cutoffs = 0.1, fdr_cutoffs = 0.4)

# Use multiple cutoffs for p and fdr
output_multistars <- strand_bias_test(strand_counts,
  p_cutoffs = c(0.5, 0.1, 0.05),
  fdr_cutoffs = c(0.5, 0.35, 0.1)
)

# Tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("tbl_df")))
  expect_true(inherits(output_repli, c("tbl_df")))
  expect_true(inherits(output_lenientcutoff, c("tbl_df")))
  expect_true(inherits(output_multistars, c("tbl_df")))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(18, 10))
  expect_equal(dim(output_repli), c(18, 10))
  expect_equal(dim(output_lenientcutoff), c(18, 10))
  expect_equal(dim(output_multistars), c(18, 10))
})

test_that("Number significant is correct", {
  expect_equal(sum(output$significant == "*"), 1)
  expect_equal(sum(output$significant_fdr == "*"), 0)
  expect_equal(sum(output_repli$significant == "*"), 0)
  expect_equal(sum(output_repli$significant_fdr == "*"), 0)
  expect_equal(sum(output_lenientcutoff$significant == "*"), 3)
  expect_equal(sum(output_lenientcutoff$significant_fdr == "*"), 3)
  expect_equal(sum(output_multistars$significant == "***"), 1)
  expect_equal(sum(output_multistars$significant_fdr == "**"), 3)
})
