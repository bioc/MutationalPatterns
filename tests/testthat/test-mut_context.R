context("test-mut_context")

#Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Get mutation context
input = unlist(vcfs)
output <- mut_context(input, ref_genome)

#Unit tests
test_that("Output has correct class",{
    expect_true(inherits(output, c("character")))
})

test_that("Output size is correct",{
    expect_equal(length(output), length(input))
})

test_that("The 64 possible contexts are returned",{
    contexts = sort(unique(output))
    expect_equal(contexts, c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", 
                               "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", 
                               "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", 
                               "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", 
                               "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", 
                               "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", 
                               "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", 
                               "TTT"))
})