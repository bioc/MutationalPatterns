context("test-binomial_test")

output_signi <- binomial_test(0.5, 1200, 543)
output_notsigni <- binomial_test(0.2, 800, 170)
output_strictcutoff <- binomial_test(0.5, 1200, 543, p_cutoffs = 0.00001)


test_that("Output has correct class", {
  expect_true(inherits(output_signi, c("data.frame")))
  expect_true(inherits(output_notsigni, c("data.frame")))
  expect_true(inherits(output_strictcutoff, c("data.frame")))
})

test_that("Output has correct size", {
  expect_equal(dim(output_signi), c(1, 3))
  expect_equal(dim(output_notsigni), c(1, 3))
  expect_equal(dim(output_strictcutoff), c(1, 3))
})

test_that("Output has correct significance level", {
  expect_equal(round(output_signi$pval, 5), 0.0011)
  expect_equal(round(output_notsigni$pval, 5), 0.39961)
  expect_equal(round(output_strictcutoff$pval, 5), 0.0011)
})

test_that("enrichment/depletion correctly determined", {
  expect_equal(output_signi$effect, factor("depletion"))
  expect_equal(output_notsigni$effect, factor("enrichment"))
  expect_equal(output_strictcutoff$effect, factor("depletion"))
})
