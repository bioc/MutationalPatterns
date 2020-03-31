#' Get contexts from larger deletions
#' 
#' @details
#' Determines the COSMIC context for deletions larger than 1bp in a GRanges object containing Indel mutations.
#' This function is called by get_indel_context_gr.
#' The function determines if there is microhomology for deletions that are not in repeat regions.
#' 
#' @param gr GRanges object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A modified version of the input gr. 
#' All variants that are not a deletion larger than 1bp are removed.
#' In each gr two columns have been added. 
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype. 
#' The subtype is the number of repeats or the microhomology length.
#' 
#' @examples 
#' 
#' @importFrom magrittr %>%
#' @family Indels
#' @seealso \code{\link{get_indel_context_gr}}
#' 
#'

get_big_dels = function(gr, mut_size, ref_genome){
    
    #Select mutations
    gr = gr[mut_size < -1]
    if (length(gr) == 0){
        return(gr)
    }
    mut_size = mut_size[mut_size < -1]
    
    #Get deleted bases
    del_bases = gr$REF %>% 
        as.character() %>% 
        substring(2)
    biggest_dels = del_bases %>%
        nchar() %>%
        max()
    flank_dist = biggest_dels * 20
    
    #Find extended sequence
    seq = get_extended_sequence(gr, flank_dist, ref_genome)
    
    #Determine nr. repeats.
    seq_z = stringr::str_replace_all(as.character(seq), del_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
    n_repeats = gsub("[^Z].*", "", as.character(seq_z)) %>% 
        nchar() + 1 #Remove all bases after the Zs and count how many bases are left. Add +1 for the deleted bases itself
    
    gr$muttype = stringr::str_c(abs(mut_size), "bp_deletion")
    gr$muttype_sub = n_repeats
    
    
    #Determine if there is microhomology for deletions that are not in repeat regions.
    pos_mh = gr$muttype_sub == 1 #There is always at least 1 'repeat', because of the deleted bases themselves.
    
    gr_repeat = gr[!pos_mh]
    gr_mh = gr[pos_mh]
    
    if (length(gr_mh) == 0){
        return(gr_repeat)
    }
    
    mut_size_mh = mut_size[pos_mh]
    del_bases_mh = del_bases[pos_mh]
    del_bases_s = strsplit(del_bases_mh, "")
    seq_s = strsplit(as.character(seq[pos_mh]), "")
    
    
    #Also check for microhomology to the left of the deletion. For this take the reverse sequence to the left of the deletion and the reverse deleted bases.
    rev_del_bases = reverse(del_bases_mh)
    rev_l_del_bases_s = strsplit(rev_del_bases, "")
    
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        l_flank = GenomicRanges::flank(gr_mh, biggest_dels)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    l_flank = l_flank %>% 
        GenomicRanges::trim() %>% 
        GenomicRanges::shift(1) #Trim the ranges that are extended beyond the actual length of the chromosome. #Add 1 base, because the first base in the granges obj is not deleted and should be used in the flank.
    rev_l_seq = BSgenome::getSeq(eval(as.symbol(ref_genome)), l_flank) %>% 
        reverse()
    rev_l_seq_s = strsplit(as.character(rev_l_seq), "")
    
    #For each mutation determine how many bases show hm
    nr_pos_mh = length(del_bases_s)
    nr_mh = vector("list", nr_pos_mh)
    for (i in 1:nr_pos_mh){
        del_bases_sample = del_bases_s[[i]]
        seq_s_sample = seq_s[[i]][1:length(del_bases_sample)]
        same = del_bases_sample == seq_s_sample
        r_nr_mh_sample = cumprod(same) %>% sum(na.rm = T) #Determine how many bases are the same before the first difference. na.rm is for when a sequence has been trimmed.
        
        l_del_bases_sample = rev_l_del_bases_s[[i]]
        l_seq_s_sample = rev_l_seq_s[[i]][1:length(l_del_bases_sample)]
        l_same = l_del_bases_sample == l_seq_s_sample
        l_nr_mh_sample = cumprod(l_same) %>% 
            sum(na.rm = T)
        
        nr_mh_sample = max(r_nr_mh_sample, l_nr_mh_sample)
        
        nr_mh[[i]] = nr_mh_sample
    }
    nr_mh = unlist(nr_mh)
    
    #Update gr when mh is indeed present
    mh_f = nr_mh > 0
    gr_mh$muttype[mh_f] = stringr::str_c(abs(mut_size_mh[mh_f]), "bp_deletion_with_microhomology")
    gr_mh$muttype_sub[mh_f] = nr_mh[mh_f]
    
    #Combine muts with and without mh
    gr = c(gr_mh, gr_repeat) %>% 
        sort()
    return(gr)
}