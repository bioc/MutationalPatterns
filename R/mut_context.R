#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @return Character vector with the context of the base substitutions
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context <- mut_context(vcfs[[1]], ref_genome)
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_context <- function(vcf, ref_genome) {
  # Check that the seqnames of the gr and ref_genome match
  .check_chroms(vcf, ref_genome)

  ranges <- GenomicRanges::resize(vcf, 3, fix = "center")

  vcf_context <- as.character(Biostrings::getSeq(
    BSgenome::getBSgenome(ref_genome),
    seqnames(vcf),
    start(vcf) - 1,
    end(vcf) + 1
  ))
  return(vcf_context)
}
