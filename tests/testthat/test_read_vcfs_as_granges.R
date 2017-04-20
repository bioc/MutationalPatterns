context("Function 'read_vcfs_as_granges'")

ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome)
library(ref_genome, character.only = TRUE)

sample_names <- c ( "colon1", "colon2", "colon3",
                   "intestine1", "intestine2", "intestine3",
                   "liver1", "liver2", "liver3" )

vcfs <- list.files (system.file("extdata", package="MutationalPatterns"),
                    pattern = ".vcf", full.names = TRUE)

test_that("loads multiple samples", {
    input <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
    expect_that(length(input), equals(9))
})

test_that("nuclear filter works", {
    input <- read_vcfs_as_granges(vcfs, sample_names, ref_genome)
    expected <- c(
        "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",
        "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14",
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22", "chrX", "chrY")

    expect_that(seqlevels(input), equals(expected))
})

test_that("autosomal filter works", {
    input <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "auto")
    expected <- c(
        "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",
        "chr8",  "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14",
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22")

    expect_that(seqlevels(input), equals(expected))
})

## test_that("sex chromosomes filter works", {
##     input <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "sex")
##     expected <- c("chrX", "chrY")

##     expect_that(seqlevels(input), equals(expected))
## })

test_that("unfiltered works", {
    # We use the reference genome that best fits the sample data here
    # to make sure the contig names automatically match.
    ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5"
    library(ref_genome, character.only = TRUE)

    input <- read_vcfs_as_granges(vcfs, sample_names, ref_genome, "none")
    expected <- seqlevels(get(ref_genome))

    proper_subset <- all(seqlevels(input) %in% expected)
    expect_equal(proper_subset, TRUE)
})
