context("test-plot_contribution")


# Load nmf data
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
  package = "MutationalPatterns"
))

## Plot the relative contribution
output <- plot_contribution(nmf_res$contribution)

## Plot the absolute contribution.
## When plotting absolute NMF results, the signatures need to be included.
output_absolute <- plot_contribution(nmf_res$contribution,
  nmf_res$signature,
  mode = "absolute"
)


## Only plot a subset of samples
output_subset <- plot_contribution(nmf_res$contribution,
  nmf_res$signature,
  mode = "absolute",
  index = c(1, 2)
)
## Flip the coordinates
output_flipcoord <- plot_contribution(nmf_res$contribution,
  nmf_res$signature,
  mode = "absolute",
  coord_flip = TRUE
)

# Use signature refitting results
fit_res <- readRDS(system.file("states/snv_refit.rds",
  package = "MutationalPatterns"
))

output_sigfit <- plot_contribution(fit_res$contribution)

## refitting results in absolute mode
output_sigfit_absolute <- plot_contribution(fit_res$contribution,
  mode = "absolute"
)

test_that("Output has correct class", {
  expect_true(inherits(output, "gg"))
  expect_true(inherits(output_absolute, "gg"))
  expect_true(inherits(output_subset, "gg"))
  expect_true(inherits(output_flipcoord, "gg"))
  expect_true(inherits(output_sigfit, "gg"))
  expect_true(inherits(output_sigfit_absolute, "gg"))
})
