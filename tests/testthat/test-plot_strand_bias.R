context("test-plot_strand_bias")

#Read stranded mut_mat
mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
                                    package="MutationalPatterns"))

## Load a reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

tissue <- c("colon", "colon", "colon",
            "intestine", "intestine", "intestine",
            "liver", "liver", "liver")

## Perform the strand bias test.
strand_counts = strand_occurrences(mut_mat_s, by=tissue)
strand_bias = strand_bias_test(strand_counts)

## Plot the strand bias.
output = plot_strand_bias(strand_bias)

#Repeat for replication bias.
mut_mat_repli <- readRDS(system.file("states/mut_mat_repli.rds",
                                     package="MutationalPatterns"))
strand_counts = strand_occurrences(mut_mat_repli, by=tissue)
strand_bias = strand_bias_test(strand_counts)
output_repli = plot_strand_bias(strand_bias)


test_that("Output has correct class",{
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_repli, c("gg")))
})