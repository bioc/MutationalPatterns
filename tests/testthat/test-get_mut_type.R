context("test-get_mut_type")

#Get a grl with variants.
grl <- readRDS(system.file("states/blood_grl.rds",
                package="MutationalPatterns"))

## Get a specific mutation type.
grl_snv = get_mut_type(grl, "snv")
grl_indel = get_mut_type(grl, "indel")
grl_dbs = get_mut_type(grl, "dbs")
grl_mbs = get_mut_type(grl, "mbs")
gr_singlesample = get_mut_type(grl[[1]], type = "dbs")
empty_gr = get_mut_type(grl[[1]][0], type = "dbs")
gr_nodbs = get_mut_type(grl[[1]][1:20], type = "dbs")

#Change names of grl_indel, to make them prettier.
remove_names_gr = function(gr){
    names(gr) = seq_along(gr)
    return(gr)
}
grl_indel = purrr::map(as.list(grl_indel), remove_names_gr) %>% 
    GRangesList()

expected_grl_indel <- readRDS(system.file("states/blood_grl_indel.rds",
                           package="MutationalPatterns"))


test_that("Output has correct class",{
    expect_true(inherits(grl_snv, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_indel, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_dbs, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(grl_mbs, c("GRanges", "CompressedGRangesList")))
    expect_true(inherits(gr_singlesample, c("GRanges")))
    expect_true(inherits(empty_gr, c("GRanges")))
    expect_true(inherits(gr_nodbs, c("GRanges")))
})

test_that("Output is equal to expected", {
    expect_equal(grl_indel, expected_grl_indel)
})

test_that("Empty gr is returned when a mut type is not present",{
    expect_equal(length(empty_gr), 0)
})

test_that("Empty gr as input results in a empty output gr",{
    expect_equal(length(gr_nodbs), 0)
})
