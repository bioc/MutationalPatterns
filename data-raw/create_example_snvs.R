library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"),
                        pattern = "sample.vcf", full.names = TRUE)

sample_names <- c(
    "colon1", "colon2", "colon3",
    "intestine1", "intestine2", "intestine3",
    "liver1", "liver2", "liver3")

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
saveRDS(grl, "inst/states/read_vcfs_as_granges_output.rds")

mut_mat = mut_matrix(grl, ref_genome)
saveRDS(mut_mat, "inst/states/mut_mat_data.rds")
