context("test-set_dbs_context")


## Get GRangesList with DBS.
grl_dbs <- readRDS(system.file("states/blood_grl_dbs.rds",
                package="MutationalPatterns"))

##Set context dbs
output = set_dbs_context(grl_dbs)

expected <- readRDS(system.file("states/blood_grl_dbs_context.rds",
                               package="MutationalPatterns"))


test_that("Output has correct class",{
    expect_true(inherits(output, c("GRanges", "CompressedGRangesList")))
})

test_that("Output is equal to expected", {
    expect_equal(output, expected)
})