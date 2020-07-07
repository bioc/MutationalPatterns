context("test-type_context")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Get type_context
input <- vcfs[[1]]
output <- type_context(input, ref_genome)
output_longer <- type_context(input, ref_genome, extension = 2)


# Unit tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("list")))
  expect_true(inherits(output$types, c("character")))
  expect_true(inherits(output$context, c("character")))
  expect_true(inherits(output_longer, c("list")))
  expect_true(inherits(output_longer$types, c("character")))
  expect_true(inherits(output_longer$context, c("character")))
})

test_that("Output size is correct", {
  expect_equal(length(output$types), length(input))
  expect_equal(length(output$context), length(input))
  expect_equal(length(output_longer$types), length(input))
  expect_equal(length(output_longer$context), length(input))
})

test_that("GRanges with 0 muts as input gives list with two empty vectors", {
  expect_warning({
    output_empty <- type_context(input[0], ref_genome)
  })
  expect_true(inherits(output_empty, "list"))
  expect_equal(length(output_empty$types), 0)
  expect_equal(length(output_empty$context), 0)
})
