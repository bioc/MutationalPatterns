setwd("hpc/hpc/cog_bioinf/cuppen/project_data/Bastiaan_signature_analysis/MutationalPatterns/R")

library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
library(NMF)
library(ggdendro)
library(pracma)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

source("read_vcfs_as_granges.R")
source("mut_96_occurrences.R")
source("mut_context.R")
source("mut_dbs_occurrences.R")
source("mut_matrix.R")
source("mut_type.R")
source("plot_profiles.R")
source("type_context.R")
source("mutations_from_vcf.R")
source("MutationalPatterns.R")
source("extract_signatures.R")
source("plot_contribution.R")
source("plot_contribution_heatmap.R")
source("plot_compare_profiles.R")
source("cos_sim.R")
source("fit_to_signatures.R")
source("cos_sim_matrix.R")
source("mut_strand.R")
source("mut_matrix_stranded.R")
source("mut_strand_occurrences.R")
source("check_mutation_type.R")
source("strand_occurrences.R")
source("strand_bias_test.R")
source("plot_strand.R")
source("plot_strand_bias.R")
source("plot_192_profile.R")
source("plot_signature_strand_bias.R")
source("binomial_test.R")
source("plot_strand_profiles.R")
source("plot_rainfall.R")

vcf_files <- c("../../../HMF_data/DR010-update/data/160704_HMFregCPCT_FR10301986_FR12244591_CPCT02020306/CPCT02020306R_CPCT02020306T_post_processed_v2.2.vcf.gz",
               "../../../HMF_data/DR010-update/data/160704_HMFregCPCT_FR12244543_FR12244602_CPCT02030245/CPCT02030245R_CPCT02030245T_post_processed_v2.2.vcf.gz",
               "../../../HMF_data/DR010-update/data/160704_HMFregCPCT_FR12244545_FR12244599_CPCT02100011/CPCT02100011R_CPCT02100011T_post_processed_v2.2.vcf.gz")

sample_names = c("t1","t2","t3")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
mode = "snv+dbs"

vcfs = read_vcfs_as_granges(vcf_files, sample_names, ref_genome, check_alleles = T)

plot_profiles(mut_matrix(vcfs, ref_genome, mode, num_cores = 1), condensed = T)

# Seperate mutation type signatures

nmf_res = extract_signatures(mut_matrix(vcfs, ref_genome, mode, num_cores = 1), rank = 2, nrun = 1)
colnames(nmf_res$signatures$snv) <- c("SNV A", "SNV B")
rownames(nmf_res$contribution$snv) <- c("SNV A", "SNV B")
colnames(nmf_res$signatures$dbs) <- c("DBS A", "DBS B")
rownames(nmf_res$contribution$dbs) <- c("DBS A", "DBS B")

plot_profiles(nmf_res$signatures, mut_type = "dbs")
plot_contribution(nmf_res$contribution, nmf_res$signatures, mode = "relative")
plot_contribution_heatmap(nmf_res$contribution$dbs, sig_order = c("SNV A", "SNV B"))
plot_compare_profiles(profile1 = mut_matrix(vcfs, ref_genome, mode, num_cores = 1)$snv, profile2 = nmf_res$reconstructed$snv)

# Combined mutation type signatures

nmf_res = extract_signatures(mut_matrix(vcfs, ref_genome, mode, method = "combine", num_cores = 1), rank = 2, nrun = 1)
colnames(nmf_res$signatures) <- c("Sig A", "Sig B")
rownames(nmf_res$contribution) <- c("Sig A", "Sig B")

plot_profiles(nmf_res$signatures)
plot_contribution(nmf_res$contribution, nmf_res$signatures, mode = "relative")
plot_contribution_heatmap(nmf_res$contribution, sig_order = c("Sig A", "Sig B"))
plot_compare_profiles(profile1 = mut_matrix(vcfs, ref_genome, mode, method = "combine", num_cores = 1)[,1], profile2 = nmf_res$reconstructed[,1])

# Find optmial contribution of COSMIC signatures

cancer_signatures = read.table("../../cosmic_signatures_v2.txt", header = T, sep = "\t")
new_order = match(TRIPLETS_96, cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures_snv = as.matrix(cancer_signatures[,4:33])

cs = read.csv("../../sigProfiler_DBS_signatures.csv")
rownames(cs) = cs$Mutation.Type
cs = as.matrix(cs[,2:11])

cancer_signatures = list("snv"=cancer_signatures_snv, "dbs"=cs)

fit_res = fit_to_signatures(mut_matrix(vcfs, ref_genome, mode, num_cores = 1), cancer_signatures)
plot_contribution(fit_res$contribution, cancer_signatures, mode = "both")
plot_contribution_heatmap(fit_res$contribution, mut_type = "snv+dbs", cluster_mut_type = "snv")
plot_compare_profiles(profile1 = mut_matrix(vcfs, ref_genome, mode, method = "split", num_cores = 1)$snv[,1], profile2 = fit_res$reconstructed$snv[,1])

cos_sim_ori_rec <- cos_sim_matrix(mut_matrix(vcfs, ref_genome, mode, method = "split", num_cores = 1), fit_res$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))

# Transcriptional strand bias

genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
strand = mut_strand(vcfs[[1]], genes_hg19, mut_type = "dbs", mode = "transcription")
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19, mut_type = "snv+dbs")
strand_counts <- strand_occurrences(mut_mat_s, mode = "snv+dbs", method = "split")
strand_bias <- strand_bias_test(strand_counts, mode = "snv+dbs")
ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand_bias(strand_bias)

plot_grid(ps1, ps2, nrow=2)

# Replication strand bias

#repli_file = system.file("../inst/extdata/ReplicationDirectionRegions.bed",package = "MutationalPatterns")
repli_file = "../inst/extdata/ReplicationDirectionRegions.bed"
repli_strand = read.table(repli_file, header = TRUE)
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1, end = repli_strand$Stop),
                               strand_info = repli_strand$Class)
seqlevelsStyle(repli_strand_granges) = "UCSC"

strand_rep <- mut_strand(vcfs[[1]], repli_strand_granges, mode = "replication")
mut_mat_s_rep <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges, mode = "replication")
strand_counts_rep <- strand_occurrences(mut_mat_s_rep)
strand_bias_rep <- strand_bias_test(strand_counts_rep)
ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)

plot_grid(ps1, ps2, nrow=2)

# Extract signatures with strand

nmf_res = extract_signatures(mut_mat_s_rep, rank = 2, nrun = 1, method = "split")
colnames(nmf_res$signatures$snv) <- c("SNV A", "SNV B")
rownames(nmf_res$contribution$snv) <- c("SNV A", "SNV B")
colnames(nmf_res$signatures$dbs) <- c("DBS A", "DBS B")
rownames(nmf_res$contribution$dbs) <- c("DBS A", "DBS B")

a <- plot_strand_profiles(nmf_res$signatures, mode = "replication", condensed = TRUE)
b <- plot_signature_strand_bias(nmf_res$signatures)

plot_grid(a,b,rel_widths = c(3,1))

# Rainfall plot

chromosomes <- seqnames(get(ref_genome))[1:22]
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), mut_type = "snv",
              chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)
