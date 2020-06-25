context("test-plot_correlation_bootstrap")

# Get contri boots
contri_boots <- readRDS(system.file("states/bootstrapped_snv_refit.rds",
  package = "MutationalPatterns"
))

# Run default function
output <- plot_correlation_bootstrap(contri_boots)

# Run for all samples combined
output_combi <- plot_correlation_bootstrap(contri_boots, per_sample = FALSE)

# Test
test_that("Output has correct class", {
  expect_true(inherits(output, c("list")))
  expect_true(inherits(output[[1]], c("gg")))
  expect_true(inherits(output_combi, c("gg")))
})
