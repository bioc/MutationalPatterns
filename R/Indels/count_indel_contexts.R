#' Count indel contexts
#' 
#' @details
#' Counts the number of indels per COSMIC context from a GRanges or GRangesList object containing Indel mutations.
#' This function applies the count_indel_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#' 
#' @param grl GRanges or GRangesList object containing Indel mutations in which the context was added with get_indel_context.
#' 
#' @return A tibble containing the number of indels per COSMIC context per gr.
#' 
#' @examples 
#' 
#' @family Indels
#' 
#' @seealso \code{\link{count_indel_contexts_gr}}, \code{\link{get_indel_context}}
#' 
#' @export


count_indel_contexts = function(grl){
    categories = tibble("muttype" = c(rep("C_deletion", 6), rep("T_deletion", 6), rep("C_insertion", 6), 
                                      rep("T_insertion", 6), rep("2bp_deletion", 6), rep("3bp_deletion", 6), 
                                      rep("4bp_deletion", 6), rep("5+bp_deletion", 6), rep("2bp_insertion", 6),
                                      rep("3bp_insertion", 6),rep("4bp_insertion", 6),rep("5+bp_insertion", 6), 
                                      rep("2bp_deletion_with_microhomology", 1), rep("3bp_deletion_with_microhomology", 2), 
                                      rep("4bp_deletion_with_microhomology", 3), rep("5+bp_deletion_with_microhomology", 5)),
                        "muttype_sub" =  c(rep(c(1:5, "6+"), 2), 
                                           rep(c(0:4, "5+"), 2), 
                                           rep(c(1:5, "6+"), 4), 
                                           rep(c(0:4, "5+"), 4), 1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        counts_l = purrr::map(gr_l, count_indel_contexts_gr, categories)
        counts = do.call(cbind, counts_l)
        colnames(counts) = names(grl)
        
    } else if (inherits(grl, "GRanges")){
        counts = count_indel_contexts_gr(grl)
        colnames(counts) = "My_sample"
    } else{
        not_gr_or_grl(grl)
    }
    counts = cbind(categories, counts)
    counts[is.na(counts)] = 0
    counts = as_tibble(counts)
    counts$muttype = factor(counts$muttype, levels = unique(counts$muttype))
    return(counts)
}