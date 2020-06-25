context("test-pool_mut_matrix")

# Get mut_mat
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))
grouping <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))


output <- pool_mut_mat(mut_mat, grouping)

test_that("Output has correct class", {
  expect_true(inherits(output, c("matrix")))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output), c(96, 3))
})
