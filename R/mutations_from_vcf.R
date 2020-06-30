#' Retrieve base substitutions from vcf
#'
#' A function to extract base substitutions of each position in vcf
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitutions
#' @import GenomicRanges
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' muts <- mutations_from_vcf(vcfs[[1]])
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mutations_from_vcf <- function(vcf) {

  # Check that no indels are present.
  .check_no_indels(vcf)

  ref <- as.character(.get_ref(vcf))
  alt <- as.character(unlist(.get_alt(vcf)))

  muts <- paste(ref, alt, sep = ">")
  return(muts)
}
