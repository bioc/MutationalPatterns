#' Removes SNV/MNV variants.
#'
#' @details 
#' This function removes SNV/MNV variants for a single GRanges object.
#' 
#' @param gr GRanges object
#' 
#' @return A filtered version of the input GRanges object.
#' @examples
#' 

remove_snvs_gr = function(gr){
    snv_f = find_snv(gr)
    gr = gr[!snv_f]
    return(gr)
}