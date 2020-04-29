context("test-plot_long_profile")

# Read the long mutation matrix information:
input <- readRDS(system.file("states/mut_mat_longregions.rds",
                                package="MutationalPatterns"))

## Plot the 96-profile of three samples
output = plot_long_profile(input)
output_condensed = plot_long_profile(input, condensed = TRUE)

test_that("Output has correct class",{
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_condensed, c("gg")))
})