context("test-binomial_test")

output_signi = binomial_test (0.5, 1200, 543)
output_notsigni = binomial_test (0.2, 800, 170)

test_that("Output has correct class",{
    expect_true(inherits(output_signi, c("data.frame")))
    expect_true(inherits(output_notsigni, c("data.frame")))
})

test_that("Output has correct size",{
    expect_equal(dim(output_signi), c(1,3))
    expect_equal(dim(output_notsigni), c(1,3))
})

test_that("Output has correct significance level",{
    expect_equal(round(output_signi$pval, 5), 0.00055)
    expect_equal(round(output_notsigni$pval, 5), 0.1998)
})

test_that("enrichment/depletion correctly determined", {
    expect_equal(output_signi$effect, factor("depletion"))
    expect_equal(output_notsigni$effect, factor("enrichment"))
})
