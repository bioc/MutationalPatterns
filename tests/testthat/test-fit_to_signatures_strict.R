context("test-fit_to_signatures_strict")

#Get mut_mat
mut_mat = readRDS(system.file("states/mut_mat_data.rds",
                              package="MutationalPatterns"))

#Get signatures
filename <- system.file("extdata/snv_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1,2)])


output = fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.05)

expected <- readRDS(system.file("states/strict_snv_refit.rds",
                                package="MutationalPatterns"))

test_that("Output has correct class",{
    expect_true(inherits(output, "list"))
    expect_true(inherits(output$fit_res, "list"))
    expect_true(inherits(output$fit_res$contribution, "matrix"))
    expect_true(inherits(output$fit_res$reconstructed, "matrix"))
    expect_true(inherits(output$sim_decay_fig, "list"))
    expect_true(inherits(output$sim_decay_fig[[1]], "gg"))
})

test_that("Output is equal to expected", {
    expect_equal(output$fit_res, expected)
})
