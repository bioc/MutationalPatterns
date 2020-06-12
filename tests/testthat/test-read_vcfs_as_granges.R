context("test-read_vcfs_as_granges")

ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(ref_genome, character.only = TRUE)

sample_names <- c ( "colon1", "colon2", "colon3",
                   "intestine1", "intestine2", "intestine3",
                   "liver1", "liver2", "liver3" )

vcfs <- list.files (system.file("extdata", package="MutationalPatterns"),
                    pattern = "sample.vcf", full.names = TRUE)

#Test default
test_that("loads multiple samples", {
    output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
    expect_that(length(output), equals(9))
    expect_true(inherits(output, "CompressedGRangesList"))
})

#Test for seqlevel filters
test_that("nuclear filter works", {
    output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
    expected <- c(
        "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",
        "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14",
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22", "chrX", "chrY")

    expect_that(seqlevels(output), equals(expected))
})

test_that("autosomal filter works", {
    output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "auto")
    expected <- c(
        "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",
        "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14",
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22")

    expect_that(seqlevels(output), equals(expected))
})

test_that("unfiltered works", {
    # We use the reference genome that best fits the sample data here
    # to make sure the contig names automatically match.
    ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5"
    library(ref_genome, character.only = TRUE)

    output <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "none")
    expected <- seqlevels(BSgenome::getBSgenome(ref_genome))

    proper_subset <- all(seqlevels(output) %in% expected)
    expect_equal(proper_subset, TRUE)
})


#Test that a warning is given when vcf and names lengths don't match
test_that("An error is given when vcf and names lengths don't match",{
    expect_error({read_vcfs_as_granges(vcfs, sample_names[1:8], ref_genome)},
                 "Please provide the same number of sample names as VCF files")
})

#Test that a warning is given when the supplied reference is not a BSgenome object
test_that("An error is given when the supplied ref is not a BSgenome object",
          {expect_error({read_vcfs_as_granges(vcfs, sample_names, "a")},
                        "Please provide the name of a BSgenome object.")
})

#Test that you can read in specific mutation types
vcf_fnames = list.files(system.file("extdata", package="MutationalPatterns"),
                        pattern = "blood.*vcf", full.names = TRUE)
sample_names = c("AC", "ACC55", "BCH")
test_that("indels work", {
    output = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "indel")
    expect_that(length(output), equals(3))
    expect_true(inherits(output, "CompressedGRangesList"))
})

test_that("dbs work", {
    output = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "dbs")
    expect_that(length(output), equals(3))
    expect_true(inherits(output, "CompressedGRangesList"))
})

    output = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "mbs")
    test_that("mbs work", {
    expect_that(length(output), equals(3))
    expect_true(inherits(output, "CompressedGRangesList"))
})

test_that("all mutation types work", {
    output = read_vcfs_as_granges(vcf_fnames, sample_names, ref_genome, type = "all")
    expect_that(length(output), equals(3))
    expect_true(inherits(output, "CompressedGRangesList"))
})

#Test function works on an empty vcf
empty_vcf = list.files(system.file("extdata", package = "MutationalPatterns"),
                       pattern = "empty.vcf", full.names = TRUE)

test_that("Empty vcf works", {
    expect_warning({output = read_vcfs_as_granges(empty_vcf, "empty", ref_genome)}, 
                   "There were 0 variants \\(before filtering\\) found in the vcf file")
    expect_true(inherits(output, "CompressedGRangesList"))
    expect_equal(length(output[[1]]), 0)
})





