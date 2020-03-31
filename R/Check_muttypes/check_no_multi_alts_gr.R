#' Checks that there are no variants with multiple alt alleles.
#'
#' @details 
#' This function checks for a single GRanges object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#' 
#' @param gr GRanges object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_multi_alts_gr = function(gr){
    alt = gr$ALT
    nr_alts = alt %>% 
        unlist() %>% 
        length()
    if (length(gr) != nr_alts){
        stop("There should not be any variants with multiple alternative alleles. Please remove these variants with remove_multi_alts_variants.", call. = F)
    }
    invisible(gr)
}