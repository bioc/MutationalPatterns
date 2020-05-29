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
#'                 package="MutationalPatterns"))
#'
#' muts = mutations_from_vcf(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mutations_from_vcf = function(vcf) {
    
    #Check that no indels are present.
    check_no_indels(vcf)
    
    # Allow both uppercase and lowercase column names.
    vcf_cols = colnames(S4Vectors::mcols(vcf))
    if ("REF" %in% vcf_cols){
        ref = as.character(vcf$REF)
    } else if ("ref" %in% vcf_cols){
        ref = as.character(vcf$ref)
    } else{
        warning("Some of your data is missing a REF column.")
        ref = character()
    }
    
    if ("ALT" %in% vcf_cols){
        alt = as.character(unlist(vcf$ALT))
    } else if ("alt" %in% vcf_cols){
        alt = as.character(unlist(vcf$alt))
    } else{
        warning("Some of your data is missing an ALT column.")
        alt = character()
    }

    muts = paste(ref, alt, sep=">")
    return(muts)
}
