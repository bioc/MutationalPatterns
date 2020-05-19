context("test-extract_signatures")

#Load mutation matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                    package="MutationalPatterns"))


#extract signatures
output <- extract_signatures(mut_mat, rank = 2, nrun = 1)

#Check it also works with indels
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
                               package="MutationalPatterns"))
output_indel = extract_signatures(indel_counts, rank = 2, nrun = 1)

#Tests
test_that("Outputs a list", {
    expect_true(inherits(output, c("list")))
    expect_true(inherits(output_indel, c("list")))
})

test_that("Output contains signatures, contribution and reconstructed", {
    expect_identical(names(output), c("signatures", "contribution", "reconstructed"))
    expect_identical(names(output_indel), c("signatures", "contribution", "reconstructed"))
})

test_that("Output elements are matrixes", {
    expect_true(inherits(output$signatures, "matrix"))
    expect_true(inherits(output$contribution, "matrix"))
    expect_true(inherits(output$reconstructed, "matrix"))
    expect_true(inherits(output_indel$signatures, "matrix"))
    expect_true(inherits(output_indel$contribution, "matrix"))
    expect_true(inherits(output_indel$reconstructed, "matrix"))
})

test_that("Output elements have correct size", {
    expect_equal(dim(output$signatures), c(96, 2))
    expect_equal(dim(output$contribution), c(2, 9))
    expect_equal(dim(output$reconstructed), c(96, 9))
    expect_equal(dim(output_indel$signatures), c(83, 2))
    expect_equal(dim(output_indel$contribution), c(2, 3))
    expect_equal(dim(output_indel$reconstructed), c(83, 3))
})