#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome. The user an choose how many bases
#' are extracted.
#'
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
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
#' ## Get the standard context
#' mut_context <- mut_context(vcfs[[1]], ref_genome)
#'
#' ## Get larger context
#' mut_context_larger <- mut_context(vcfs[[1]], ref_genome, extension = 2)
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_context <- function(vcf, ref_genome, extension = 1) {
  # Check that the seqnames of the gr and ref_genome match
  .check_chroms(vcf, ref_genome)

  # Get context of mutation.
  vcf_context <- as.character(Biostrings::getSeq(
    BSgenome::getBSgenome(ref_genome),
    seqnames(vcf),
    start(vcf) - extension,
    end(vcf) + extension
  ))
  return(vcf_context)
}
