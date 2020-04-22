context("test-fit_to_signatures")

#Get mut_mat
mut_mat = readRDS(system.file("states/mut_mat_data.rds",
                              package="MutationalPatterns"))

#Get signatures
filename <- system.file("extdata/snv_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1,2)])


output = fit_to_signatures(mut_mat, signatures)

expected <- readRDS(system.file("states/snv_refit.rds",
                                package="MutationalPatterns"))

test_that("Output has correct class",{
    expect_true(inherits(output, "list"))
    expect_true(inherits(output$contribution, "matrix"))
    expect_true(inherits(output$reconstructed, "matrix"))

})

test_that("Output is equal to expected", {
    expect_equal(output, expected)
})
