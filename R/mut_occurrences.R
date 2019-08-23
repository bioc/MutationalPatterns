#' Count double base substitutions occurrences
#'  
#' @param type_context result of double base substitutions from type_context function
#' @importFrom S4Vectors isEmpty
#' @noRd
#' @return vector with 78 double base substitutions mutation occurrences

mut_occurrences = function(type_context, mode, indel)
{
  mode = check_mutation_type(mode)
  
  count = list()
 
  for (m in mode)
  {
    # if type_context is empty, return vector with zeroes for mutation type
    if (all(isEmpty(type_context[[m]])))
    {
      count[[m]] = NULL
      next
    }
    
    if (m == "snv")
    {
      count[[m]] = rep(0,96)
      names(count[[m]]) = TRIPLETS_96
      context = sprintf("%s[%s]%s", 
                        substr(type_context[[m]][["context"]], 1, 1),
                        type_context[[m]][["types"]],
                        substr(type_context[[m]][["context"]], 3, 3))
      context = table(context)
    } else if (m == "dbs") 
    {
      count[[m]] = rep(0,78)
      names(count[[m]]) = DBS
      context = table(type_context[[m]][["types"]])
    } else if (m == "indel")
    {
      if (missing(indel)) { indel = "native" }
      if (indel == "cosmic")
      { 
        count[[m]] = rep(0,83) 
        names(count[[m]]) = indel_context
      } else if (indel == "native")
      { 
        count[[m]]=rep(0,30)
        names(count[[m]]) = indel_context
      }
      
      context = table(type_context[[m]][["context"]])
    }
    
    count[[m]][names(context)] = context
  }
  
  return(count)
}

##
## Deprecated variants
##

#'
#' This function has been removed.  Use 'mut_occurrences' instead.
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

mut_96_occurrences <- function(vcf, ref_genome)
{
  .Defunct("mut_occurrences", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_occurrences' instead."))
}
