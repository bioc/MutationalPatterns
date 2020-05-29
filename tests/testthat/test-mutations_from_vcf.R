context("test-mutations_from_vcf")

#Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

output = mutations_from_vcf(vcfs[[1]])
output_empty = mutations_from_vcf(vcfs[[1]][0])


#Unit tests
test_that("Output has correct class",{
    expect_true(inherits(output, c("character")))
    expect_true(inherits(output_empty, c("character")))
})

test_that("The 12 substitution types are returned",{
    types = sort(unique(output))
    expect_equal(types, c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", 
                               "G>A", "G>C", "G>T","T>A", "T>C", "T>G"))
})

test_that("GRanges with 0 muts as input gives empty output",{
    expect_equal(length(output_empty), 0)
})
