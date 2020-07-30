#' Retrieve context of base substitution types
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitution types.
#'
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @return Mutation types and context character vectors in a named list
#'
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
#' ## Get type context
#' type_context <- type_context(vcfs[[1]], ref_genome)
#'
#' ## Get larger type context
#' type_context_larger <- type_context(vcfs[[1]], ref_genome, extension = 2)
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_context}}
#'
#' @export

type_context <- function(vcf, ref_genome, extension = 1) {
  # Deal with empty GRanges objects.
  if (length(vcf) == 0) {
    warning("Detected empty GRanges object.
                Returning an empty list for this sample.", call. = FALSE)
    res <- list("types" = NULL, "context" = NULL)
    return(res)
  }

  # Get the mut context
  mut_context <- mut_context(vcf, ref_genome, extension)

  # Get the mutations
  muts <- mutations_from_vcf(vcf)

  # Get the 6 base mutation types
  types <- mut_type(vcf)

  # find the mutations for which the context needs to be adjusted
  x <- which(muts != types)

  # subset mut_context
  y <- mut_context[x]

  # Change the context of these mutations to reverse complement
  # of the context
  y <- IRanges::reverse(chartr("ATGC", "TACG", y))

  # replace subset with reverse complement
  mut_context[x] <- y

  # return as named list
  res <- list(types, mut_context)
  names(res) <- c("types", "context")

  return(res)
}
