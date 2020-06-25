context("test-mutations_from_vcf")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))
vcf <- vcfs[[1]]

# Run function
output <- mutations_from_vcf(vcf)

# Check it works on empty input
output_empty <- mutations_from_vcf(vcf[0])

# Check it works on lowercase input
vcf_lowercase <- vcf
colnames(mcols(vcf_lowercase)) <- c("paramRangeID", "ref", "alt", "QUAL", "FILTER")
output_lowercase <- mutations_from_vcf(vcf_lowercase)

# Check it gives a warning on data with no ref or alt
vcf_noref <- vcf
colnames(mcols(vcf_noref)) <- c("paramRangeID", "a", "ALT", "QUAL", "FILTER")
vcf_noalt <- vcf
colnames(mcols(vcf_noalt)) <- c("paramRangeID", "REF", "a", "QUAL", "FILTER")

# Unit tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("character")))
  expect_true(inherits(output_empty, c("character")))
})

test_that("The 12 substitution types are returned", {
  types <- sort(unique(output))
  expect_equal(types, c(
    "A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"
  ))
})

test_that("GRanges with 0 muts as input gives empty output", {
  expect_equal(length(output_empty), 0)
})

test_that("Input with lowercase doesn't change result", {
  expect_equal(output, output_lowercase)
})

test_that("GR with no REF or ALT gives an error.", {
  expect_error(
    {
      output_noref <- mutations_from_vcf(vcf_noref)
    },
    "missing a REF column"
  )
  expect_error(
    {
      output_noalt <- mutations_from_vcf(vcf_noalt)
    },
    "missing a ALT column"
  )
})
