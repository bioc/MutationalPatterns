#' Count the occurrences of each base substitution type
#'
#' @param vcf_list GRangesList or GRanges object.
#' @param ref_genome BSGenome reference genome object
#' @return data.frame with counts of each base substitution type for
#' each sample in vcf_list
#'
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_type_occurrences <- function(vcf_list, ref_genome) {


  # Convert to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  } else if (inherits(vcf_list, "GRanges")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
    names(vcf_list) <- "my_sample"
  }

  # Check that the seqnames of the gr and ref_genome match
  .check_chroms(vcf_list, ref_genome)

  # Check input
  if (!inherits(vcf_list, "CompressedGRangesList")) {
    .not_gr_or_grl(vcf_list)
  }

  n_samples <- length(vcf_list)
  df <- data.frame()

  CpG <- c("ACG", "CCG", "TCG", "GCG")
  column_names <- c(
    "C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
    "C>T at CpG", "C>T other"
  )

  full_table <- NULL
  for (i in seq_len(n_samples))
  {
    vcf <- vcf_list[[i]]
    types <- mut_type(vcf)

    CT_context <- 0
    CT_at_CpG <- 0
    CT_at_other <- 0

    CT_muts <- which(types == "C>T")
    if (length(CT_muts) > 0) {
      CT_context <- type_context(vcf[CT_muts], ref_genome)[[2]]
      CT_at_CpG <- sum(!(is.na(BiocGenerics::match(CT_context, CpG))))
      CT_at_other <- length(CT_muts) - CT_at_CpG
    }

    # Construct a table and handle missing mutation types.
    full_table <- table(factor(types, levels = column_names))
    full_table["C>T at CpG"] <- CT_at_CpG
    full_table["C>T other"] <- CT_at_other
    df <- BiocGenerics::rbind(df, full_table)
  }

  row.names(df) <- names(vcf_list)
  colnames(df) <- names(full_table)
  return(df)
}
