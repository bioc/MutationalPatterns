# Use biomaRt to obtain data.

mart <- "ensemble"
library(biomaRt)
regulatory <- useEnsembl(
  biomart = "regulation",
  dataset = "hsapiens_regulatory_feature",
  GRCh = 37
)
saveRDS(regulatory, "inst/states/regulatory_data.rds")

# Download the regulatory CTCF binding sites and convert them to
# a GRanges object.
CTCF <- getBM(
  attributes = c(
    "chromosome_name",
    "chromosome_start",
    "chromosome_end",
    "feature_type_name"
  ),
  filters = "regulatory_feature_type_name",
  values = "CTCF Binding Site",
  mart = regulatory
)

CTCF_g <- reduce(GRanges(
  CTCF$chromosome_name,
  IRanges(
    CTCF$chromosome_start,
    CTCF$chromosome_end
  )
))
seqlevels(CTCF_g) <- c(1:22, "X", "Y")
CTCF_g <- sort(CTCF_g)
CTCF_g <- CTCF_g[sample.int(length(CTCF_g), 50000)]
saveRDS(CTCF_g, "inst/states/CTCF_g_data.rds")

# Download the promoter regions and conver them to a GRanges object.
promoter <- getBM(
  attributes = c(
    "chromosome_name", "chromosome_start",
    "chromosome_end", "feature_type_name"
  ),
  filters = "regulatory_feature_type_name",
  values = "Promoter",
  mart = regulatory
)
promoter_g <- reduce(GRanges(
  promoter$chromosome_name,
  IRanges(
    promoter$chromosome_start,
    promoter$chromosome_end
  )
))
seqlevels(promoter_g) <- c(1:22, "X", "Y")
promoter_g <- sort(promoter_g)
saveRDS(promoter_g, "inst/states/promoter_g_data.rds")

flanking <- getBM(
  attributes = c(
    "chromosome_name",
    "chromosome_start",
    "chromosome_end",
    "feature_type_name"
  ),
  filters = "regulatory_feature_type_name",
  values = "Promoter Flanking Region",
  mart = regulatory
)
flanking_g <- reduce(GRanges(
  flanking$chromosome_name,
  IRanges(
    flanking$chromosome_start,
    flanking$chromosome_end
  )
))
seqlevels(flanking_g) <- c(1:22, "X", "Y")
flanking_g <- sort(flanking_g)
flanking_g <- flanking_g[sample.int(length(flanking_g), 50000)]
saveRDS(flanking_g, "inst/states/promoter_flanking_g_data.rds")
