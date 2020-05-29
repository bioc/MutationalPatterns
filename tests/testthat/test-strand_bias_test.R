context("test-strand_bias_test")

#Load stranded mutation matrix
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
                                    package="MutationalPatterns"))

## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Set tissue names
tissue <- c("colon", "colon", "colon",
            "intestine", "intestine", "intestine",
            "liver", "liver", "liver")

## Perform the strand bias test.
strand_counts = strand_occurrences(mut_mat_s, by=tissue)
output = strand_bias_test(strand_counts)

#Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
                                     package="MutationalPatterns"))
strand_counts = strand_occurrences(mut_mat_repli, by=tissue)
output_repli = strand_bias_test(strand_counts)


#Tests
test_that("Output has correct class",{
    expect_true(inherits(output, c("tbl_df")))
    expect_true(inherits(output_repli, c("tbl_df")))
})

test_that("Output has correct size",{
    expect_equal(dim(output), c(18, 8))
    expect_equal(dim(output_repli), c(18, 8))
})