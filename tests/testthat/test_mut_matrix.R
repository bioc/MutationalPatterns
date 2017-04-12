context("Function 'mut_matrix'")

# To test mut_matrix, we need to load the reference genome first.
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# We re-use the data that is shipped with the package.
sample_names <- c ( "colon1", "colon2", "colon3",
                   "intestine1", "intestine2", "intestine3",
                   "liver1", "liver2", "liver3" )

vcfs <- list.files (system.file("extdata", package="MutationalPatterns"),
                    pattern = ".vcf", full.names = TRUE)

input <- read_vcfs_as_granges(vcfs, sample_names,
                              "BSgenome.Hsapiens.UCSC.hg19")

expected <- readRDS("mut_matrix.rds")

test_that("transforms correctly", {
    actual <- mut_matrix(input, ref_genome)
    expect_that(actual, equals(expected))
})

test_that("a list and a GRangesList are acceptable", {
    list_actual <- mut_matrix(as.list(input), ref_genome)
    grangeslist_actual <- mut_matrix(input, ref_genome)

    expect_that (list_actual, equals(grangeslist_actual))
    expect_that (grangeslist_actual, equals(expected))
})
