library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


#Read vcfs into granges
vcf_fnames = list.files("inst/extdata/", pattern = "blood.*.vcf", full.names = T)

vcf_l = purrr::map(vcf_fnames, readVcf, "hg19")
sample_names = basename(vcf_fnames) %>% 
    str_remove("blood-") %>% 
    str_remove(".vcf")
names(vcf_l) = sample_names
grl = purrr::map(vcf_l, granges) %>% 
    GRangesList()
seqlevelsStyle(grl) = "UCSC"
seqlevels(grl, pruning.mode = "fine") = str_c("chr", c(1:22, "X"))
saveRDS(grl, "inst/states/blood_grl.rds")

#Get indels
grl_indel = get_mut_type(grl, "indel")
saveRDS(grl_indel, "inst/states/blood_grl_indel.rds")

grl_indel_context = get_indel_context(grl_indel, ref_genome)
saveRDS(grl_indel_context, "inst/states/blood_grl_indel_context.rds")

indel_counts = count_indel_contexts(grl_indel_context)
saveRDS(indel_counts, "inst/states/blood_indel_counts.rds")


