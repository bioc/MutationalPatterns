context("test-cos_sim")

# Calculate cosine similarity
x <- c(1.1, 2.1, 0.2, 0.1, 2.9)
y <- c(0.9, 1.9, 0.5, 0.4, 3.1)
output <- cos_sim(x, y)


test_that("Output has correct class and data type", {
  expect_true(inherits(output, c("numeric")))
  expect_equal(typeof(output), "double")
})

test_that("Output has expected size", {
  expect_equal(length(output), 1)
})

test_that("Output is equal to expected", {
  expect_equal(output, 0.9895599)
})
