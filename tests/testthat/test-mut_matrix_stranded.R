context("test-mut_matrix_stranded")


# To test mut_matrix, we need to load the reference genome and the genes first.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

# Test that the function works with default arguments
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
input <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))
expected <- readRDS(system.file("states/mut_mat_s_data.rds",
  package = "MutationalPatterns"
))

test_that("transforms correctly", {
  output <- mut_matrix_stranded(input, ref_genome, ranges = genes_hg19)
  expect_equal(output, expected)
})

# Test that a list is an acceptable input
test_that("a list is also acceptable input", {
  output <- mut_matrix_stranded(input, ref_genome, ranges = genes_hg19)
  output_list <- mut_matrix_stranded(as.list(input), ref_genome, ranges = genes_hg19)

  expect_equal(output_list, output)
  expect_equal(output_list, expected)
})

# A single sample can be used as input.
test_that("A single GR can also be used as input", {
  output_singlesample <- mut_matrix_stranded(input[[1]], ref_genome, ranges = genes_hg19)
  expect_true(inherits(output_singlesample, "matrix"))
  expect_equal(dim(output_singlesample), c(192, 1))
})

# seqlevels genes need to match the input
genes_badseqlevel <- genes_hg19
seqlevels(genes_badseqlevel)[1] <- "chrtest"
test_that("A single GR can also be used as input", {
  expect_error(
    {
      mut_matrix_stranded(input[[1]], ref_genome, ranges = genes_badseqlevel)
    },
    "Chromosome names \\(seqlevels\\) of vcf and genes Granges object do not match"
  )
})


# Test replication mode
repli_strand_granges <- readRDS(system.file("states/repli_strand.rds",
  package = "MutationalPatterns"
))
expected_repli <- readRDS(system.file("states/mut_mat_repli.rds",
  package = "MutationalPatterns"
))

test_that("replication mode transforms correctly", {
  mut_mat_repli <- mut_matrix_stranded(input, ref_genome, repli_strand_granges, mode = "replication")
  expect_equal(mut_mat_repli, expected_repli)
})


# Test longer context
output_longer <- mut_matrix_stranded(input, ref_genome, ranges = genes_hg19, extension = 2)

test_that("Output has correct class", {
  expect_true(inherits(output_longer, "matrix"))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output_longer), c(3072, 9))
})
