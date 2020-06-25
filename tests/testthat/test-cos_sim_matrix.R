context("test-cos_sim_matrix")

# Read signatures
filename <- system.file("extdata/snv_signatures_probabilities.txt",
  package = "MutationalPatterns"
)
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures <- as.matrix(signatures[, -c(1, 2)])

# Read mut_matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))


# Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
output <- cos_sim_matrix(mut_mat, signatures)

test_that("Output has correct class and data type", {
  expect_true(inherits(output, c("matrix")))
  expect_equal(typeof(output), "double")
})

test_that("Output has expected size", {
  expect_equal(dim(output), c(9, 47))
})
