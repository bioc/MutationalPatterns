context("test-plot_compare_profiles")

#Load mutation data
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                                package="MutationalPatterns"))

#Load nmf data
nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
                    package="MutationalPatterns"))

#Compare profiles
output = plot_compare_profiles(mut_mat[,1],
                        nmf_res$reconstructed[,1],
                        profile_names = c("Original", "Reconstructed"))

output_condensed = plot_compare_profiles(mut_mat[,1],
                               nmf_res$reconstructed[,1],
                               profile_names = c("Original", "Reconstructed"),
                               condensed = TRUE)

#Perform tests
test_that("Output has correct class",{
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_condensed, c("gg")))
})
