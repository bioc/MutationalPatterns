#' Checks that there are no SNVs/MNVs.
#'
#' @details 
#' This function checks for a single GRanges object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#' 
#' @param gr GRanges object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_snvs_gr = function(gr){
    snv_f = find_snv(gr)
    nr_snv = sum(snv_f)
    snv_present = nr_snv >= 1
    if (snv_present){
        stop(str_c("There seem to be ", nr_snv, " SNVs present in your data. Make sure to remove all SNVs with the remove_snvs function."), call. = F)
    }
    invisible(gr)
}