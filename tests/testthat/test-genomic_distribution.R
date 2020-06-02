context("test-genomic_distribution")

#Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                package="MutationalPatterns"))

#Load genomic regions
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
                    package="MutationalPatterns"))

promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
                        package="MutationalPatterns"))

flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
                                    package="MutationalPatterns"))

#Combine regions and set seqlevelstyle
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
seqlevelsStyle(regions) <- "UCSC"

#Get the callable regions
surveyed_file <- system.file("extdata/callableloci-sample.bed",
                            package="MutationalPatterns")

library(rtracklayer)
surveyed <- import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"

#Use the same callable loci for all samples.
surveyed_list <- rep(list(surveyed), 9)

## Calculate the number of observed and expected number of mutations in
## each genomic regions for each sample.
output <- genomic_distribution(vcfs, surveyed_list, regions)


test_that("Output has correct class",{
    expect_true(inherits(output, c("data.frame")))
})

test_that("Output has correct size",{
    expect_equal(dim(output), c(27, 8))
})
