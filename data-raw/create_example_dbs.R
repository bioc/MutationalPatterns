library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Get grl
grl = readRDS("inst/states/blood_grl.rds")

#Get indels
grl_dbs = get_mut_type(grl, "dbs")
saveRDS(grl_dbs, "inst/states/blood_grl_dbs.rds")

grl_dbs = set_context_dbs(grl_dbs)
saveRDS(grl_dbs, "inst/states/blood_grl_dbs_context.rds")
