library(tidyverse)
library(VariantAnnotation)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Get grl
grl <- readRDS("inst/states/blood_grl.rds")

# Get mbs
grl_mbs <- get_mut_type(grl, "mbs")
saveRDS(grl_mbs, "inst/states/blood_grl_mbs.rds")

# Count contexts
mbs_counts <- count_mbs_contexts(grl_mbs)
saveRDS(mbs_counts, "inst/states/blood_mbs_counts.rds")
