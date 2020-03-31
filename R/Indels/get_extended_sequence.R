#' Retreive the flanking sequence from a GRanges object
#' 
#' @details
#' This function retreives the flanking sequence from a GRanges object.
#' The size of the flanking sequence is supplied as an argument to the function.
#' This function works by first extending the ranges and then retreiving the sequence.
#' 
#'
#' 
#' @param gr GRanges object containing Indel mutations. 
#' The mutations should be called similarly to HaplotypeCaller.
#' @param flank_dist A numeric vector of length one containing the number of flanking base pairs.
#' @param ref_genome BSGenome reference genome object
#' 
#' @return A DNAStringSet containing the flaning bases.
#' 
#' @examples 
#' 
#' @family Indels
#'
#' 
#'

get_extended_sequence = function(gr, flank_dist, ref_genome){
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_extended = GenomicRanges::flank(gr, flank_dist, start = F)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_extended = GenomicRanges::trim(gr_extended) #Trim the ranges that are extended beyond the actual length of the chromosome.
    seq = BSgenome::getSeq(eval(as.symbol(ref_genome)), gr_extended)
    return(seq)
}