context("test-plot_profile_heatmap")


# Get input data
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

mut_mat_extended <- readRDS(system.file("states/mut_mat_data_extended.rds",
  package = "MutationalPatterns"
))


## Create heatmap of profile
output_basic <- plot_profile_heatmap(mut_mat, max = 0.1)

## Create heatmap of extended profile
output <- plot_profile_heatmap(mut_mat_extended)

## Or plot heatmap per tissue
tissue <- c(
  "colon", "colon", "colon",
  "intestine", "intestine", "intestine",
  "liver", "liver", "liver"
)

output_tissue <- plot_profile_heatmap(mut_mat_extended, by = tissue)

## Or plot the heatmap per sample.
output_sample <- plot_profile_heatmap(mut_mat_extended,
  by = colnames(mut_mat_extended),
  max = 0.05
)


test_that("Output has correct class", {
  expect_true(inherits(output_basic, c("gg")))
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_tissue, c("gg")))
  expect_true(inherits(output_sample, c("gg")))
})
