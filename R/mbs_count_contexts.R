#' Count MBS variants grouped by length.
#' 
#' @details
#' Counts the number of mbs grouped by length from a GRanges or GRangesList object containing mbs variants.
#' This is used, since a COSMIC context has to our knowledge not yet been defined.
#' This function applies the count_mbs_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#' 
#' @param grl GRanges or GRangesList object containing mbs variants.
#' 
#' @return A tibble containing the number of MBS per MBS length per gr.
#' 
#' @examples 
#' ## Get a GRangesList or GRanges object with mbs variants.
#' grl_mbs <- readRDS(system.file("states/blood_grl_mbs.rds",
#'                 package="MutationalPatterns"))
#' 
#' #Count the MBSs
#' count_mbs_contexts(grl_mbs)
#' 
#' @family MBS
#' @seealso \code{\link{count_mbs_contexts_gr}}
#' 
#' @export
count_mbs_contexts = function(grl){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    size = NULL
    
    categories = tibble::tibble("size" = c(3:9, "10+"))
    
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        counts_l = purrr::map(gr_l, count_mbs_contexts_gr, categories)
        counts = do.call(cbind, counts_l)
        colnames(counts) = names(grl)
        
    } else if (inherits(grl, "GRanges")){
        counts = count_mbs_contexts_gr(grl, categories)
        colnames(counts) = "My_sample"
    } else{
        not_gr_or_grl(grl)
    }
    counts = cbind(categories, counts)
    
    #Turn output into matrix
    counts = counts %>% 
        tibble::column_to_rownames("size") %>% 
        as.matrix()
    return(counts)
}

#' Count MBS grouped by length from a single GRanges object.
#' 
#' @details
#' Counts the number of MBS per MBS length from a GRanges object containing mbs variants.
#' The function is called by count_mbs_contexts
#' 
#' @param gr GRanges object containing mbs variants.
#' @param categories A tibble containing all possible mbs size categories
#' 
#' @return A single column tibble containing the number of MBS per MBS length
#' 
#' @importFrom magrittr %>%
#' 
#' @seealso \code{\link{count_mbs_contexts}}
#' @family MBS
count_mbs_contexts_gr = function(gr, categories){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    count = size = . = NULL
    
    #Create count table
    counts_tb = gr %>% #Determine different sizes
        get_ref() %>% 
        BiocGenerics::width() %>% 
        tibble::enframe(value = "size") %>% 
        dplyr::select(-name) %>%
        dplyr::mutate(size = ifelse(size >= 10, "10+", size),
                      size = as.character(size)) %>% 
        dplyr::group_by(size) %>% #Count number muts per size
        dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>% 
        dplyr::left_join(categories, ., by = "size") %>% #Add to possible categories
        dplyr::mutate(count = ifelse(is.na(count), 0, count)) %>% 
        dplyr::select(-size)
    
    return(counts_tb)
}
