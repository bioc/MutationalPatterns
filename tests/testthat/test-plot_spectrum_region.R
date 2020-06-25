context("test-plot_spectrum_region")

# load data
grl <- readRDS(system.file("states/grl_split_region.rds",
  package = "MutationalPatterns"
))

# Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


# Get the type occurrences for all VCF objects.
type_occurrences <- mut_type_occurrences(grl, ref_genome)

# Plot the point mutation spectrum over all samples
output <- plot_spectrum_region(type_occurrences)

# Plot the point mutation spectrum, relative only to the samples.
output_relative_sample <- plot_spectrum_region(type_occurrences, mode = "relative_sample")


# Plot the absolute point mutation spectrum over all samples
output_absolute <- plot_spectrum_region(type_occurrences, mode = "absolute")

# Plot using different types of error bars.
output_stdev <- plot_spectrum_region(type_occurrences, error_bars = "stdev")
output_sem <- plot_spectrum_region(type_occurrences, error_bars = "SEM")

# Plot including the sample points themselves
output_indv <- plot_spectrum_region(type_occurrences, indv_points = TRUE)

# Plot per tissue
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)
output_tissue <- plot_spectrum_region(type_occurrences, by = tissue)

# Plot each sample separately
sample_names <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)
output_sample <- plot_spectrum_region(type_occurrences, by = sample_names, error_bars = "none")

# Test different outputs
test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_relative_sample, c("gg")))
  expect_true(inherits(output_absolute, c("gg")))
  expect_true(inherits(output_stdev, c("gg")))
  expect_true(inherits(output_sem, c("gg")))
  expect_true(inherits(output_indv, c("gg")))
  expect_true(inherits(output_tissue, c("gg")))
  expect_true(inherits(output_sample, c("gg")))
})
