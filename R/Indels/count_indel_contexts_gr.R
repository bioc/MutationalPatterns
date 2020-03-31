#' Count indel contexts from a single gr.
#' 
#' @details
#' Counts the number of indels per COSMIC context from a GRanges object containing Indel mutations.
#' The function is called by count_indel_contexts_gr
#' 
#' @param gr GRanges object containing Indel mutations in which the context was added with get_indel_context.
#' @param categories A tibble containing all possible indel context categories
#' 
#' @return A single column tibble containing the number of indels per COSMIC context.
#' 
#' @examples 
#' 
#' @importFrom magrittr %>%
#' @family Indels
#' 
#' @seealso \code{\link{count_indel_contexts}}, \code{\link{get_indel_context}}


count_indel_contexts_gr = function(gr, categories){
    if (length(gr) == 0){
        categories = categories %>%
            dplyr::mutate(count = 0) %>% 
            dplyr::select(-muttype, -muttype_sub)
        return(categories)
    }
    
    id_context = tibble("muttype" = gr$muttype, "muttype_sub" = gr$muttype_sub)
    
    #Classify the number of repeat units/ homopolymer length / microhomology length to either 5+ or 6+ depending on whether the indel is a ins or del.
    id_context[id_context$muttype_sub >= 6, "muttype_sub"] = "6+"
    id_context[grepl("insertion|microhomology", id_context$muttype) & id_context$muttype_sub >= 5, "muttype_sub"] = "5+"
    
    #Classify large indels as size 5+
    ref_sizes = gr$REF %>%
        width()
    alt_sizes = gr$ALT %>% 
        unlist() %>% 
        width()
    mut_size = abs(alt_sizes - ref_sizes)
    mut_size_f = mut_size >= 5
    id_context$muttype = ifelse(mut_size_f, gsub("[0-9]+bp", "5+bp", id_context$muttype, perl = T), id_context$muttype)
    
    id_context_count = id_context %>% 
        dplyr::group_by(muttype, muttype_sub) %>% 
        dplyr::summarise(count = dplyr::n())
    id_context_count_full = dplyr::left_join(categories, id_context_count, by = c("muttype", "muttype_sub")) %>% 
        dplyr::select(-muttype, -muttype_sub)
    #colnames(id_context_count_full)[3] = name
    return(id_context_count_full)
}