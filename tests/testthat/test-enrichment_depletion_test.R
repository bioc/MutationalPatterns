context("test-enrichment_depletion_test")

#Read distribution data
distr <- readRDS(system.file("states/distr_data.rds",
                    package="MutationalPatterns"))
#Set tissue
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))

## Perform the enrichment/depletion test by tissue type.
output <- enrichment_depletion_test(distr, by = tissue)

## Or without specifying the 'by' parameter.
output_singlesample <- enrichment_depletion_test(distr)

test_that("Output has correct class",{
    expect_true(inherits(output, c("data.frame")))
    expect_true(inherits(output_singlesample, c("data.frame")))
})

test_that("Output has correct size",{
    expect_equal(dim(output), c(15, 13))
    expect_equal(dim(output_singlesample), c(45,13))
})
