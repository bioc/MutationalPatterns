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

expected <- matrix(c(6, 7, 4, 4, 8, 3, 2, 3, 9, 3,
                     1, 2, 7, 5, 2, 11, 4, 3, 0, 4,
                     1, 2, 0, 1, 1, 0, 0, 0, 3, 1,
                     0, 3, 14, 7, 64, 10, 5, 12, 29, 15,
                     10, 12, 38, 3, 10, 6, 16, 12, 2, 4,
                     4, 4, 0, 0, 3, 2, 0, 0, 0, 0,
                     1, 0, 2, 4, 13, 2, 3, 16, 2, 3,
                     2, 6, 4, 2, 2, 3, 4, 3, 3, 9,
                     0, 0, 4, 3, 1, 2, 0, 2, 1, 2,
                     0, 1, 2, 1, 1, 5,6, 3, 2, 3,
                     4, 2, 1, 6, 5, 4, 4, 5, 6, 3,
                     2, 12, 3, 1, 0, 1, 1, 2, 0, 2,
                     0, 0, 0, 0, 3, 3, 0, 1, 11, 5,
                     63, 6, 13, 5, 29, 15, 8, 14, 27, 4,
                     8, 13, 25, 13, 8, 3, 2, 5, 2, 3,
                     3, 1, 0, 0, 0, 1, 1, 1, 3, 5,
                     9, 2, 6, 12, 1, 7, 5, 7, 3, 3,
                     4, 7, 4, 3, 3, 8, 1, 1, 2, 4,
                     1, 0, 1, 1, 1, 2, 0, 0, 2, 1,
                     2, 3,5, 9, 1, 4, 9, 2, 0, 0,
                     10, 3, 0, 1, 2, 7, 2, 13, 2, 2,
                     1, 1, 0, 2, 1, 1, 1, 0, 0, 2,
                     0, 4, 0, 5, 10, 5, 53, 10, 9, 6,
                     34, 10, 7, 12, 45, 7, 8, 11, 24, 10,
                     2, 3, 6, 5, 2, 1, 5, 2, 1, 0,
                     0, 2, 1, 1, 1, 12, 14, 4, 2, 10,
                     2, 4, 3, 4, 6, 5, 2, 2, 4, 1,
                     6, 4, 2, 2, 0, 1, 0, 0, 2, 1,
                     0, 0, 1, 1, 1, 1, 4,3, 10, 7,
                     2, 2, 7, 3, 1, 8, 8, 1, 2, 3,
                     8, 6, 3, 8, 2, 1, 3, 2, 0, 1,
                     0, 3, 2, 0, 1, 1, 2, 3, 0, 4,
                     14, 2, 60, 6, 7, 9, 23, 13, 11, 8,
                     42, 5, 8, 10, 22, 11, 4, 3, 1, 4,
                     0, 0, 1, 2, 1, 0, 1, 0, 2, 3,
                     0, 6, 10, 5, 6, 11, 2, 3, 8, 7,
                     2, 2, 7, 2, 5, 3, 5, 9, 1, 1,
                     1, 0, 3, 0, 2, 0, 0, 1, 0, 2,
                     3, 1, 1, 2, 10, 7, 1, 5, 6, 4,
                     1, 3, 10, 4, 3, 3, 8, 3, 1, 13,
                     4, 1, 0, 1, 2, 2, 1, 3, 2, 0,
                     0, 0, 2, 0, 0, 6, 15, 6, 48, 8,
                     9, 12, 31, 9, 9, 7, 32, 8, 6, 14,
                     28, 14, 1, 5, 7, 4, 3, 0, 0, 1,
                     0, 0, 0, 0, 3, 4, 1, 5, 13, 3,
                     6, 9, 1, 3, 3, 5, 4, 1, 3, 4,
                     4, 3, 4, 4, 0, 0, 5, 1, 1, 1,
                     1, 1, 0, 1, 0, 0, 1, 1, 2, 3,
                     10, 8, 0, 2, 9, 5, 0, 3, 5, 2,
                     1, 6, 7, 3, 3, 8, 3, 1, 2, 2,
                     1, 2, 0, 2, 4, 0, 0, 1, 4, 3,
                     0, 3, 11, 9, 58, 11, 8, 11, 29, 10,
                     7, 6, 43, 7, 7, 7, 23, 5, 4, 3,
                     4, 2, 1, 2, 2, 2, 1, 0, 1, 1,
                     2, 3, 1, 6, 8, 3, 3, 12, 4, 6,
                     5, 8, 6, 3, 2, 3, 4, 2, 3, 5,
                     2, 0, 1, 4, 2, 0, 4, 0, 1, 2,
                     0, 0, 2, 0, 0, 4, 12, 7, 0, 4,
                     9, 5, 0, 1, 8, 3, 3, 4, 8, 5,
                     1, 15, 4, 0, 1, 2, 3, 2, 2, 3,
                     0, 0, 1, 1, 6, 1, 1, 5, 14, 8,
                     48, 8, 9, 7, 22, 12, 7, 11, 41, 11,
                     13, 7, 21, 12, 4, 3, 1, 3, 2, 0, 4,
                     1, 1, 0, 1, 2, 2, 3, 0, 5, 5,
                     2, 5, 10, 4, 6, 4, 4, 4, 2, 2,
                     6, 3, 3, 1, 10, 0, 1, 2, 2, 0,
                     0, 2, 0, 0, 0, 0, 0, 2, 1, 2,
                     5, 9, 6, 3, 2, 4, 2, 1, 1, 9,
                     2, 3, 6, 5, 2, 2, 6, 1, 1, 1,
                     3, 0, 1, 1, 2, 1, 0, 0, 0, 1,
                     1, 0, 5, 13, 9, 60, 6, 10, 5, 28,
                     13, 14, 10, 31, 11, 10, 12, 27, 10, 4,
                     0, 4, 7, 1, 1, 3, 3, 1, 0, 0,
                     0, 1, 1, 1, 7, 11, 2, 8, 13, 3,
                     1, 3, 5, 7, 3, 3, 6, 9, 1, 5,
                     2, 0, 1, 2, 0, 0, 1, 1, 1, 0,
                     1, 3, 2, 4, 3, 1, 3,  9, 3, 2,
                     3, 10, 2, 1, 5, 7, 3, 5, 4, 6,
                     3, 2, 9, 2, 3, 0, 4, 0, 3, 0,
                     3, 1, 1, 1, 0, 1, 2, 0, 6, 8,
                     5, 68, 6, 10, 7, 32, 15, 4, 10, 33,
                     12, 12, 4, 22, 12, 4, 4, 3, 2, 0,
                     0, 0, 2, 0, 1, 1, 2, 0, 2, 4,
                     5, 10, 0, 7, 7, 4, 4, 10, 10, 4,
                     0, 2, 5, 4, 2, 4, 7, 1, 1, 1,
                     2, 1, 0, 1, 2, 1, 0, 1, 0, 3,
                     1, 2, 2),
                   nrow=96,
                   ncol=9)

dimnames(expected) <- list(c("ACA", "ACC", "ACG", "ACT",
                             "CCA", "CCC", "CCG", "CCT",
                             "GCA", "GCC", "GCG", "GCT",
                             "TCA", "TCC", "TCG", "TCT",
                             "ACA", "ACC", "ACG", "ACT",
                             "CCA", "CCC", "CCG", "CCT",
                             "GCA", "GCC", "GCG", "GCT",
                             "TCA", "TCC", "TCG", "TCT",
                             "ACA", "ACC", "ACG", "ACT",
                             "CCA", "CCC", "CCG", "CCT",
                             "GCA", "GCC", "GCG", "GCT",
                             "TCA", "TCC", "TCG", "TCT",
                             "ATA", "ATC", "ATG", "ATT",
                             "CTA", "CTC", "CTG", "CTT",
                             "GTA", "GTC", "GTG", "GTT",
                             "TTA", "TTC", "TTG", "TTT",
                             "ATA", "ATC", "ATG", "ATT",
                             "CTA", "CTC", "CTG", "CTT",
                             "GTA", "GTC", "GTG", "GTT",
                             "TTA", "TTC", "TTG", "TTT",
                             "ATA", "ATC", "ATG", "ATT",
                             "CTA", "CTC", "CTG", "CTT",
                             "GTA", "GTC", "GTG", "GTT",
                             "TTA", "TTC", "TTG", "TTT"),
                           c("colon1", "colon2", "colon3",
                             "intestine1", "intestine2", "intestine3",
                             "liver1", "liver2", "liver3"))


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
