#' Retrieve context of mutations
#' 
#' A function to extract the contexts of mutations for all mutation types
#'
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @param mode A character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' @return Mutation types and context character vectors in a named list
#'
#' @importFrom IRanges reverse
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
#' type_context <- type_context(vcfs[[1]], ref_genome, mode)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_context}}
#'
#' @export

type_context = function(vcf, ref_genome, mode, ...)
{
    mode = check_mutation_type(mode)
    
    # Deal with empty GRanges objects.
    if (length (vcf) == 0)
    {
        warning("Detected empty GRanges object.")
        res = list("types"=c(), "context"=c())
        return(res)
    }
    
    res = list()

    for (m in mode)
    {
      muts = mutations_from_vcf(vcf, m)
      types = mut_type(vcf, m)
      
      # if snvs are extracted, convert the base substitutions to the
      # conventional base substitution types
      if (m == "snv")
      {
        mut_context = mut_context(vcf, ref_genome, m)$snv
        
        # find the mutations for which the context needs to be adjusted
        x = which(muts$snv != types$snv)
        
        # subset mut_context
        y = mut_context[x]
        
        # Change the context of these mutations to reverse complement
        # of the context
        y = reverse(chartr('ATGC', 'TACG', y))
        
        # replace subset with reverse complement
        mut_context[x] = y
        
        res[[m]] = list("types"=types[[m]], "context"=mut_context)
      } else if (m == "dbs"){
        res[[m]] = list("types"=types[[m]])
      } else if (m == "indel"){
        mut_context = mut_context(vcf, ref_genome, m, ...)$indel
        res[[m]] = list("types"=types[[m]], "context"=mut_context)
      }
    }
    
    if (length(res) == 1)
      res = list("context"=res[[1]][["context"]], "types"=res[[1]][["types"]])
    # return as named list
    return(res)
}