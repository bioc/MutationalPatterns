#' Retrieve base substitution types from a VCF object
#' 
#' A function to extract the base substitutions from a vcf and translate to
#' the 6 common base substitution types.
#' 
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitution types
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' mutation_types(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mutation_types = function(vcf) 
{
    muts = mutations_from_vcf(vcf)
    types = unlist(muts)
    types = gsub('G>T', 'C>A', types)
    types = gsub('G>C', 'C>G', types)
    types = gsub('G>A', 'C>T', types)
    types = gsub('A>T', 'T>A', types)
    types = gsub('A>G', 'T>C', types)
    types = gsub('A>C', 'T>G', types)
    return(types)
}

##
## Deprecated variants
##

get_types <- function(vcf)
{
    .Defunct("mutation_types", package="MutationalPatterns",
                msg=paste("This function has been removed.  Use",
                            "'mutation_types' instead."))
}
