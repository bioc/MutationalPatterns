#' Removes variants with multiple alt alleles.
#'
#' @details 
#' This function removes variants with multiple alternative alleles for a single GRanges object.
#' 
#' @param gr GRanges object
#' 
#' @return A filtered version of the input GRanges object.
#' @examples
#' 

remove_multi_alts_variants_gr = function(gr){
    alt = gr$ALT
    gr = gr[elementNROWS(alt) == 1]
    return(gr)
}