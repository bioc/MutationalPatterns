#' Get indel contexts
#' 
#' @details
#' Determines the COSMIC context from a GRanges or GRangesList object containing Indel mutations.
#' It applies the get_indel_context_gr function to each gr in the input.
#' 
#' @param grl GRanges or GRangesList object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A modified version of the input grl. In each gr two columns have been added. 
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype. 
#' The subtype is either the number of repeats or the microhomology length. 
#' 
#' @examples 
#' 
#' @family Indels
#' 
#' @seealso
#' \code{\link{get_indel_context_gr}},
#' 
#' @export

get_indel_context = function(grl, ref_genome){
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        gr_list = purrr::map(gr_l, function(x) get_indel_context_gr(x, ref_genome))
        grl = GenomicRanges::GRangesList(gr_list)
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = get_indel_context_gr(grl, ref_genome)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}