library(tidyverse)
library(VariantAnnotation)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Get grl
grl = readRDS("inst/states/blood_grl.rds")

#Get indels
grl_indel = get_mut_type(grl, "indel")
saveRDS(grl_indel, "inst/states/blood_grl_indel.rds")

#Get context
grl_indel_context = get_indel_context(grl_indel, ref_genome)
saveRDS(grl_indel_context, "inst/states/blood_grl_indel_context.rds")

#Count contexts
indel_counts = count_indel_contexts(grl_indel_context)
saveRDS(indel_counts, "inst/states/blood_indel_counts.rds")


#Refit to signatures
filename <- system.file("extdata/indel_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1)])

indel_m = indel_counts %>% 
    dplyr::select(-muttype, -muttype_sub) %>% 
    as.matrix()
fit_res = fit_to_signatures(indel_m, signatures)
saveRDS(fit_res, "inst/states/indel_refit.rds")

