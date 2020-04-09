context("test-get_mut_type")

#Get a grl with variants.
grl <- readRDS(system.file("states/blood_grl.rds",
                package="MutationalPatterns"))

## Get a specific mutation type.
grl_snv = get_mut_type(grl, "snv")
grl_indel = get_mut_type(grl, "indel")
grl_dbs = get_mut_type(grl, "dbs")
grl_mbs = get_mut_type(grl, "mbs")

expected_grl_indel <- readRDS(system.file("states/blood_grl_indel.rds",
                           package="MutationalPatterns"))


test_that("Output has correct class",{
    expect_true(inherits(grl_snv, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_indel, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_dbs, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_mbs, c("GRanges", "CompressedGRangesList")))
})

test_that("Output is equal to expected", {
    expect_equal(grl_indel, expected_grl_indel)
})
