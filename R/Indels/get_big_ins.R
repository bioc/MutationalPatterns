#' Get contexts from bigger inserions
#' 
#' @details
#' Determines the COSMIC context for insertions larger than 1bp in a GRanges object containing Indel mutations.
#' This function is called by get_indel_context_gr.
#'
#' 
#' @param gr GRanges object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A modified version of the input gr. 
#' All variants that are not insertions larger than 1bp are removed.
#' In each gr two columns have been added. 
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype. 
#' The subtype is the number of repeats.
#' 
#' @examples 
#' 
#' @family Indels
#' @seealso \code{\link{get_indel_context_gr}}
#' 
#'

get_big_ins = function(gr, mut_size, ref_genome){
    
    #Select mutations
    gr = gr[mut_size > 1]
    if (length(gr) == 0){
        return(gr)
    }
    mut_size = mut_size[mut_size > 1]
    
    #Get inserted bases
    ins_bases = gr$ALT %>% 
        unlist() %>% 
        as.character() %>% 
        substring(2)
    biggest_ins = ins_bases %>% 
        nchar() %>% 
        max()
    flank_dist = biggest_ins * 20
    
    #Get extended sequence
    seq = get_extended_sequence(gr, flank_dist, ref_genome)
    
    #Determine nr. repeats.
    seq_z = str_replace_all(as.character(seq), ins_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
    n_repeats = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() #Remove all bases after the Zs and count how many bases are left.
    
    #Return results
    gr$muttype = str_c(mut_size, "bp_insertion")
    gr$muttype_sub = n_repeats
    return(gr)
}