#' Get indel contexts from a single gr
#' 
#' @details
#' Determines the COSMIC context from a GRanges object containing Indel mutations.
#' It throws an error if there are any variants with multiple alternative alleles or SNVs.
#' 
#' @param gr GRanges object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A modified version of the input gr. In the gr two columns have been added. 
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype. 
#' The subtype is either the number of repeats or the microhomology length. 
#' 
#' @examples 
#' 
#' @family Indels
#' 
#' @seealso \code{\link{get_indel_context}},


get_indel_context_gr = function(gr, ref_genome){
    
    #Check that no snvs are present.    
    check_no_snvs(gr)
    
    #Calculate indel size to determine main category
    ref_sizes = gr$REF %>%
        width()
    alt_sizes = gr$ALT %>%
        unlist() %>%
        width()
    mut_size = alt_sizes - ref_sizes
    
    #For the main indel categories, determine their sub categories. (Also split the big deletion categorie into repeat and micro homology.)
    gr_1b_dels = get_1bp_dels(gr, mut_size, ref_genome)
    gr_1b_ins = get_1bp_ins(gr, mut_size, ref_genome)
    gr_big_dels = get_big_dels(gr, mut_size, ref_genome)
    gr_big_ins = get_big_ins(gr, mut_size, ref_genome)
    
    gr = c(gr_1b_dels, gr_1b_ins, gr_big_dels, gr_big_ins) %>%
        sort()
    return(gr)
}