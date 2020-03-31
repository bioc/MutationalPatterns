#' Checks that there are no SNVs/MNVs.
#'
#' @details 
#' This function checks for a GRanges/GRangesList object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#' 
#' @param grl GRanges/GrangesList object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 

check_no_snvs = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr = unlist(grl)
    } else if (inherits(grl, "GRanges")){
        gr = grl
    } else{
        not_gr_or_grl(grl)
    }
    check_no_snvs_gr(gr)
    invisible(grl)
}