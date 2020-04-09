library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Determine vcf and donor names
snv_vcf_fnames = list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/hg19/Healthy_bone_marrow/", 
           pattern = "MQ60.vcf", 
           full.names = T)
indel_vcf_fnames = list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/INDELs/hg19/Healthy_bone_marrow/", 
                        pattern = "CALLABLE.vcf", 
                        full.names = T)
vcf_fnames = c(snv_vcf_fnames, indel_vcf_fnames)
sample_names = vcf_fnames %>% 
    basename() %>%
    str_remove("_.*")
donor_names = sample_names %>% 
    str_remove("HSC.*") %>% 
    str_remove("MPP.*")

used_donors = c("AC", "ACC55","BCH")
vcf_f = donor_names %in% used_donors
vcf_fnames = vcf_fnames[vcf_f]
donor_names = donor_names[vcf_f]

#Read vcfs and combine vcfs from the same donor
vcf_l = purrr::map(vcf_fnames, readVcf, genome = "hg19")
vcf_l = vcf_l %>% 
    split(donor_names) %>% #Split per donor into list of lists
    purrr::map(., function(x) do.call(rbind, x)) %>% #Combine vcfs from a single donor
    purrr::map(sort) %>% 
    purrr::map(unique)

#Write vcfs
out_vcf_fnames = str_c("inst/extdata/blood-", names(vcf_l), ".vcf")
purrr::map2(vcf_l, out_vcf_fnames, writeVcf)


#Create granges object from the vcfs
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
