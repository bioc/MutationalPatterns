context("test-set_context_dbs")


## Get GRangesList with DBS.
grl_dbs <- readRDS(system.file("states/blood_grl_dbs.rds",
                package="MutationalPatterns"))

##Set context dbs
output = set_context_dbs(grl_dbs)

expected <- readRDS(system.file("states/blood_grl_dbs_context.rds",
                               package="MutationalPatterns"))


test_that("Output has correct class",{
    expect_true(inherits(output, c("GRanges", "CompressedGRangesList")))
})

test_that("Output is equal to expected", {
    expect_equal(output, expected)
})