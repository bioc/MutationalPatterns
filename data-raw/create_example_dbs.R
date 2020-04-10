library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Get grl
grl = readRDS("inst/states/blood_grl.rds")

#Get indels
grl_dbs = get_mut_type(grl, "dbs")
saveRDS(grl_dbs, "inst/states/blood_grl_dbs.rds")

#Set context
grl_dbs_context = set_context_dbs(grl_dbs)
saveRDS(grl_dbs_context, "inst/states/blood_grl_dbs_context.rds")

#Count contexts
dbs_counts = count_dbs_contexts(grl_dbs_context)
saveRDS(dbs_counts, "inst/states/blood_dbs_counts.rds")
