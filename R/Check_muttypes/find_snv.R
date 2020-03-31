#' Identifies SNVs/MNVs
#'
#' @details 
#' This function finds SNVs/MNVs for a single GRanges object.
#' 
#' 
#' @param gr GRanges object
#' 
#' @return A boolean vector. It's TRUE for SNVs/MNVs and FALSE for Indels
#' @examples
#' 

find_snv = function(gr){
    check_no_multi_alts(gr)
    snv_f = width(gr$REF) == width(unlist(gr$ALT))
    return(snv_f)
}