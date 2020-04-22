library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

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

genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
mut_mat_s = mut_matrix_stranded(grl, ref_genome, genes_hg19)
saveRDS(mut_mat_s, "inst/states/mut_mat_s_data.rds")

repli_file = system.file("extdata/ReplicationDirectionRegions.bed",
                         package = "MutationalPatterns")
repli_strand = read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1,
                                                end = repli_strand$Stop),
                               strand_info = repli_strand$Class)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) = "UCSC"
saveRDS(repli_strand_granges, "inst/states/repli_strand.rds")

mut_mat_repli = mut_matrix_stranded(grl, ref_genome, repli_strand_granges, mode = "replication")
saveRDS(mut_mat_repli, "inst/states/mut_mat_repli.rds")
