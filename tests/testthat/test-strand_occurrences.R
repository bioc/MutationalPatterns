context("test-strand_occurrences")

# Read in stranded mutation matrix
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Set tissue names
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

output <- strand_occurrences(mut_mat_s, by = tissue)

# Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
  package = "MutationalPatterns"
))
output_repli <- strand_occurrences(mut_mat_repli, by = tissue)

# Tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("tbl_df")))
  expect_true(inherits(output_repli, c("tbl_df")))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(36, 5))
  expect_equal(dim(output_repli), c(36, 5))
})
