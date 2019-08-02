#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @param mode A character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' @return Character vector with the context of the base substitutions
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings getSeq
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context <- mut_context(vcfs[[1]], ref_genome, mode)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_context = function(vcf, ref_genome, mode="snv") 
{
    # Make sure that the chromosome names are compatible with each other.
    if (!(all(seqlevels(vcf) %in% seqlevels(get(ref_genome)))))
        stop(paste( "The chromosome names (seqlevels) of the VCF and the",
                    "reference genome object do not match. Use the",
                    "'seqlevelsStyle()' function to rename chromosome",
                    "names.") )

    if (!grepl("snv",tolower(mode)) & tolower(mode) != "all")
        stop(paste( "Extract context is only available for single base substitutions"))
  
    vcf <- vcf[nchar(as.character(vcf$REF))==1 & nchar(as.character(unlist(vcf$ALT)))==1]
  
    vcf_context = as.character(getSeq(get(ref_genome),
                                        seqnames(vcf),
                                        start(vcf) - 1,
                                        end(vcf) + 1))
    return(vcf_context)
}

##
## Deprecated variants
##

#'
#' This function has been removed.  Use 'mut_context' instead.
#'
#' @param vcf        A GRanges object
#' @param ref_genome The reference genome
#'
#' @return Character vector with the context of the base substitutions
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context <- mut_context(vcfs[[1]], ref_genome)
#'
#' @seealso
#' \code{\link{mut_context}}
#'
#' @export

mutation_context <- function(vcf, ref_genome)
{
  .Defunct("mut_context", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_context' instead."))
}
