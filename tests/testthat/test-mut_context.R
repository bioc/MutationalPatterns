context("test-mut_context")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

## Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Get mutation context
input <- unlist(vcfs)
output <- mut_context(input, ref_genome)
output_long <- mut_context(input, ref_genome, extension = 2)

# Unit tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("character")))
  expect_true(inherits(output_long, c("character")))
})

test_that("Output size is correct", {
  expect_equal(length(output), length(input))
  expect_equal(length(output_long), length(input))
})

test_that("The 64 possible contexts are returned", {
  contexts <- sort(unique(output))
  expect_equal(contexts, c(
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA",
    "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC",
    "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG",
    "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT",
    "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA",
    "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC",
    "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG",
    "TTT"
  ))
})


# Test that the wrong genome name gives an error
input_wronggenomename <- input
genome(input_wronggenomename) <- "hg38"
test_that("Wrong genome name gives error", {
  expect_error(
    {
      mut_context(input_wronggenomename, ref_genome)
    },
    "The input GRanges \\(your vcf data\\) and the ref_genome do not have the same genome name"
  )
})

# Test that the wrong seqlevelstyle gives an error.
input_wrongseqstyle <- input
genome(input_wrongseqstyle) = NA
seqlevelsStyle(input_wrongseqstyle) <- "NCBI"
genome(input_wrongseqstyle) = "hg19"
test_that("Wrong seqlevelsStyle gives error", {
  expect_error(
    {
      mut_context(input_wrongseqstyle, ref_genome)
    },
    "The input GRanges and the ref_genome share no seqnames"
  )
})

# Test that seqlevels not present in the ref gives an error
input_wrongseqlevel <- input
seqlevels(input_wrongseqlevel) <- paste0("chr", c(1:22, "X", "Y", "test"))
genome(input_wrongseqlevel) <- "hg19"
test_that("Wrong seqlevel gives error", {
  expect_error(
    {
      mut_context(input_wrongseqlevel, ref_genome)
    },
    paste0(
      "seqlevels \\(chromosome names\\) occur in the input GRanges,",
      "but are not present in the ref_genome"
    )
  )
})


# Test that variant that doesn't overlap with the reference give an error
input_nooverlap <- input
suppressWarnings({
  ranges(input_nooverlap)[1] <- IRanges(start = 1000000000, end = 1000000000)
})
test_that("Non overlapping variant gives error", {
  expect_error(
    {
      mut_context(input_nooverlap, ref_genome)
    },
    paste0(
      "variants that don't overlap with ",
      "the chromosome lengths of the chosen reference genome"
    )
  )
})
