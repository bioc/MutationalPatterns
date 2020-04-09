context("test-indel_get_context")

## Get a GRangesList object with only indels.
grl_indel <- readRDS(system.file("states/blood_grl_indel.rds",
                package="MutationalPatterns"))

## Load the corresponding reference genome.
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

## Get the indel contexts
output = get_indel_context(grl_indel, ref_genome)

expected = readRDS(system.file("states/blood_grl_indel_context.rds",
                               package="MutationalPatterns"))

test_that("Output has correct class",{
    expect_true(inherits(output, c("GRanges", "CompressedGRangesList")))
})

test_that("Output is equal to expected", {
    expect_equal(output, expected)
})
