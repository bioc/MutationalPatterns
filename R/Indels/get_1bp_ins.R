#' Get contexts from 1bp insertions
#' 
#' @details
#' Determines the COSMIC context for the 1bp insertions in a GRanges object containing Indel mutations.
#' This function is called by get_indel_context_gr.
#'
#' 
#' @param gr GRanges object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A modified version of the input gr. 
#' All variants that are not a 1bp insertion are removed.
#' In each gr two columns have been added. 
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype. 
#' The subtype is the number of repeats.
#' 
#' @examples 
#' 
#' @importFrom magrittr %>%
#' @family Indels
#' @seealso \code{\link{get_indel_context_gr}}
#' 
#'

get_1bp_ins = function(gr, mut_size, ref_genome){
    
    #Select mutations
    gr = gr[mut_size == 1]
    if (length(gr) == 0){
        return(gr)
    }
    
    #Get inserted bases.
    ins_bases = gr$ALT %>% 
        unlist() %>% 
        as.character() %>% 
        substring(2)
    
    #Get extended sequence
    seq = get_extended_sequence(gr, 20, ref_genome)
    
    #Check homopolymer length
    seq_z = stringr::str_replace_all(as.character(seq), ins_bases, rep("Z", length(seq))) #For each mut replace the inserted basetype in the flanking sequence with Zs.
    homopolymer_length = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() #Remove all bases after the Zs and count how many bases are left.
    ins_bases[ins_bases == "A"] = "T"
    ins_bases[ins_bases == "G"] = "C"
    
    #Return the results
    gr$muttype = stringr::str_c(ins_bases, "_insertion")
    gr$muttype_sub = homopolymer_length
    return(gr)
}