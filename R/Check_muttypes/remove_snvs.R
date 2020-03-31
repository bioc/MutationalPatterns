#' Removes SNV/MNV variants.
#'
#' @details 
#' This function removes SNV/MNV variants for a GRanges/GrangesList object.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return A filtered version of the input GRanges/GrangesList object.
#' @examples
#' 
#' @export

remove_snvs = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        grl = purrr::map(gr_l, remove_snvs_gr) %>% 
            GRangesList()
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = remove_snvs_gr(grl)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}