context("test-plot_strand")

# Load stranded mutation matrix data
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

# Calculate strand counts
strand_counts <- strand_occurrences(mut_mat_s, by = tissue)

# Plot the strand in relative mode.
output <- plot_strand(strand_counts)

# Plot in absolute mode.
output_absolute <- plot_strand(strand_counts, mode = "absolute")

# Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
  package = "MutationalPatterns"
))
strand_counts <- strand_occurrences(mut_mat_repli, by = tissue)
output_repli <- plot_strand(strand_counts)


test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_absolute, c("gg")))
  expect_true(inherits(output_repli, c("gg")))
})
