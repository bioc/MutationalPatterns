context("test-plot_lesion_segregation")


# Load GRangesList
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                           package="MutationalPatterns"))

#Select sample
gr = grl[[1]]

#Perform function
output = plot_lesion_segregation(gr, sample_name = "Colon1")
output_noname = plot_lesion_segregation(gr)
output_l = plot_lesion_segregation(gr, per_chrom = TRUE, sample_name = "Colon1")

test_that("Output has correct class",{
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_noname, c("gg")))
    expect_true(inherits(output_l, c("list")))
    expect_true(inherits(output_l[[1]], c("gg")))
})

test_that("Output per chromosome has correct length",{
    expect_equal(length(output_l), 22)
})