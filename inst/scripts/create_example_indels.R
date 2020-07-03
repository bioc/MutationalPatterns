library(tidyverse)
library(VariantAnnotation)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Get grl
grl <- readRDS("inst/states/blood_grl.rds")

# Get indels
grl_indel <- get_mut_type(grl, "indel")

# Remove names from gr, because they are often very long.
remove_names_gr <- function(gr) {
  names(gr) <- seq_along(gr)
  return(gr)
}
grl_indel <- purrr::map(as.list(grl_indel), remove_names_gr) %>%
  GRangesList()

saveRDS(grl_indel, "inst/states/blood_grl_indel.rds")

# Get context
grl_indel_context <- get_indel_context(grl_indel, ref_genome)
saveRDS(grl_indel_context, "inst/states/blood_grl_indel_context.rds")

# Count contexts
indel_counts <- count_indel_contexts(grl_indel_context)
saveRDS(indel_counts, "inst/states/blood_indel_counts.rds")


# Refit to signatures
signatures <- get_known_signatures("indel")

fit_res <- fit_to_signatures(indel_counts, signatures)
saveRDS(fit_res, "inst/states/indel_refit.rds")




# Split per region
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
  package = "MutationalPatterns"
))
promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
  package = "MutationalPatterns"
))
flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
  package = "MutationalPatterns"
))

# Combine the regions into a single GRangesList
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")

seqlevelsStyle(regions) <- "UCSC"
grl_indel_split <- split_muts_region(grl_indel_context, regions)
indel_counts_split <- count_indel_contexts(grl_indel_split)
saveRDS(indel_counts_split, "inst/states/blood_indels_counts_split_region.rds")
indel_matrix_long <- lengthen_mut_matrix(indel_counts_split)
saveRDS(indel_matrix_long, "inst/states/blood_indels_longmatrix_split_region.rds")
