#' Checks that there are no variants with multiple alt alleles.
#'
#' @details 
#' This function checks for a GRanges/GRangesList object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_multi_alts = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr = unlist(grl)
    } else if (inherits(grl, "GRanges")){
        gr = grl
    } else{
        not_gr_or_grl(grl)
    }
    check_no_multi_alts_gr(gr)
    invisible(grl)
}