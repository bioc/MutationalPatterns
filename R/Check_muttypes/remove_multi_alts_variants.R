#' Removes variants with multiple alt alleles.
#'
#' @details 
#' This function removes variants with multiple alternative alleles for a GRanges/GRangesList object.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return A filtered version of the input GRanges/GRangesList object.
#' @examples
#' 
#' @export

remove_multi_alts_variants = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        grl = purrr::map(gr_l, remove_multi_alts_variants_gr) %>% 
            GRangesList()
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = remove_multi_alts_variants_gr(grl)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}