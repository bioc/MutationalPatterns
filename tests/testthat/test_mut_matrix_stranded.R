context("Function 'mut_matrix_stranded'")

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"),
                        pattern = ".vcf", full.names = TRUE)

sample_names <- c("colon1", "colon2", "colon3",
                  "intestine1", "intestine2", "intestine3",
                  "liver1", "liver2", "liver3")

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

expected <- readRDS(system.file("states/mut_mat_s_data.rds",
                                package="MutationalPatterns"))

test_that("transforms correctly", {
    actual <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
    expect_that(actual, equals(expected))
})


test_that("a list and a GRangesList are acceptable", {
    list_actual <- mut_matrix_stranded(as.list(vcfs), ref_genome, genes_hg19)
    grangeslist_actual <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
    
    expect_that(list_actual, equals(grangeslist_actual))
    expect_that(grangeslist_actual, equals(expected))
})

