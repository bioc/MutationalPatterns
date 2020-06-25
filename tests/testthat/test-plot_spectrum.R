context("test-plot_spectrum")

# Laad variants:
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))


## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

## Get the type occurrences for all VCF objects.
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

## Plot the point mutation spectrum over all samples
output <- plot_spectrum(type_occurrences)

## CT distinction
output_CT <- plot_spectrum(type_occurrences, CT = TRUE)

## You can also include individual sample points.
output_indv <- plot_spectrum(type_occurrences, CT = TRUE, indv_points = TRUE)

## You can also change the type of error bars
output_stdev <- plot_spectrum(type_occurrences, error_bars = "stdev")
output_sem <- plot_spectrum(type_occurrences, error_bars = "SEM")

## Or plot spectrum per tissue
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

output_tissue <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE)

## Or plot the spectrum per sample. Error bars are set to 'none', because they can't be plotted.
output_sample <- plot_spectrum(type_occurrences, by = names(vcfs), CT = TRUE, error_bars = "none")

## You can also set custom colors.
my_colors <- c(
  "pink", "orange", "blue", "lightblue",
  "green", "red", "purple"
)

## And use them in a plot.
output_color <- plot_spectrum(type_occurrences,
  CT = TRUE,
  legend = TRUE,
  colors = my_colors
)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_CT, c("gg")))
  expect_true(inherits(output_indv, c("gg")))
  expect_true(inherits(output_stdev, c("gg")))
  expect_true(inherits(output_sem, c("gg")))
  expect_true(inherits(output_tissue, c("gg")))
  expect_true(inherits(output_sample, c("gg")))
  expect_true(inherits(output_color, c("gg")))
})
