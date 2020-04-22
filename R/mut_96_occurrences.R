#' Count 96 trinucleotide mutation occurrences
#'  
#'  @details 
#'  This function is called by mut_matrix. It calculates the 96 trinucleotide context for all variants
#'  and then splits these per GRanges (samples). It then calculates how often each 96 trinucleotide context occures.
#'  
#'  
#' @param type_context result from type_context function
#' @param gr_sizes A vector indicating the number of variants per GRanges
#' @return Mutation matrix with 96 trinucleotide mutation occurrences
#' 
#' @importFrom magrittr %>% 
mut_96_occurrences = function(type_context, gr_sizes){
    
    #Create table with all possible contexts
    cats = tibble::tibble("categories" = factor(TRIPLETS_96, levels = TRIPLETS_96))
    
    #Determine 96 context for all variants
    full_context = stringr::str_c(substr(type_context$context, 1, 1), 
                                  "[", 
                                  type_context$types, 
                                  "]", 
                                  substr(type_context$context, 3, 3)) %>% 
        factor(levels = TRIPLETS_96)
    
    #Set names if they are not yet present
    if (is.null(names(gr_sizes))){
        names(gr_sizes) = seq_along(gr_sizes)
    }
    
    #Create vector describing the sample of each variant
    sample_vector = rep(names(gr_sizes), gr_sizes) %>% 
        factor(levels = names(gr_sizes))
    
    #Count the mutations per type and per sample
    counts = tibble::tibble("categories" = full_context, "sample" = sample_vector) %>% 
        dplyr::filter(!is.na(categories)) %>% 
        dplyr::group_by(categories, sample, .drop = F) %>% 
        dplyr::summarise(count = dplyr::n())

    #Transform the data into a mutation matrix
    counts = tidyr::spread(counts, key = sample, value = count, fill = 0)
    unnecesary_cols = which(colnames(counts) == "<NA>")
    mut_mat = as.matrix(counts[,-c(1, unnecesary_cols)])
    rownames(mut_mat) = counts$categories
    return(mut_mat)
}