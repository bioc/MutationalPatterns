context("test-convert_sigs_to_ref")

# Load mutation matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

# Get signatures
signatures <- get_known_signatures(source = "SIGNAL", sig_type = "tissue", tissue_type = "Skin")

# Fit tissue specific signatures
fit_res <- fit_to_signatures(mut_mat, signatures)

# Convert the tissue specific signatures exposures to reference
output <- convert_sigs_to_ref(fit_res)

# Run tests
test_that("Output has correct class", {
  expect_true(inherits(output, "list"))
  expect_true(inherits(output$contribution, "matrix"))
  expect_true(inherits(output$reconstructed, "matrix"))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output$contribution), c(38, 9))
})

test_that("Nr. mutations hasn't changed", {
  expect_equal(colSums(output$contribution), colSums(fit_res$contribution))
})

# Test that an error is thrown when the sig names don't match.
fit_res_badname <- fit_res
rownames(fit_res_badname$contribution)[1] <- "fakename"
test_that("An error is thrown when the sig names don't match", {
  expect_error(
    {
      convert_sigs_to_ref(fit_res_badname)
    },
    "The signature names of the fit_res don't match that of"
  )
})
