#' Count DBS contexts
#' 
#' @details
#' Counts the number of dbs per COSMIC context from a GRanges or GRangesList object containing dbs variants.
#' This function applies the count_dbs_contexts_gr function to each gr in its input.
#' It then combines the results in a single tibble and returns this.
#' 
#' @param grl GRanges or GRangesList object containing DBS mutations in which the context was added with set_dbs_context.
#' 
#' @return A tibble containing the number of DBS per COSMIC context per gr.
#' 
#' @examples 
#' ## Get a GRangesList or GRanges object with dbs contexts.
#' ## See 'dbs_get_context' for more info on how to do this.
#' grl_dbs_context <- readRDS(system.file("states/blood_grl_dbs_context.rds",
#'                 package="MutationalPatterns"))
#' 
#' #Count the dbs contexts
#' count_dbs_contexts(grl_dbs_context)
#' 
#' @family DBS
#' @seealso \code{\link{count_dbs_contexts_gr}}, \code{\link{set_dbs_context}}
#' 
#' @export
count_dbs_contexts = function(grl){
    categories = tibble::tibble("REF" = c(rep("AC", 9),rep("AT", 6), rep("CC", 9), rep("CG", 6), 
                                  rep("CT", 9), rep("GC", 6), rep("TA", 6), rep("TC", 9), 
                                  rep("TG", 9), rep("TT", 9)), 
                        "ALT" = c("CA", "CG", "CT", "GA", "GG", "GT", "TA", "TG", "TT",
                                  "CA", "CC", "CG", "GA", "GC", "TA", "AA", "AG", "AT", 
                                  "GA", "GG", "GT", "TA", "TG", "TT", "AT", "GC", "GT", 
                                  "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG", 
                                  "TA", "TC", "TG", "AA", "AG", "AT", "CA", "CG", "TA", 
                                  "AT", "CG", "CT", "GC", "GG", "GT", "AA", "AG", "AT", 
                                  "CA", "CG", "CT", "GA", "GG", "GT", "AA", "AC", "AT", 
                                  "CA", "CC", "CT", "GA", "GC", "GT", "AA", "AC", "AG", 
                                  "CA", "CC", "CG", "GA", "GC", "GG"))
    
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        counts_l = purrr::map(gr_l, count_dbs_contexts_gr, categories)
        counts = do.call(cbind, counts_l)
        colnames(counts) = names(grl)
        
    } else if (inherits(grl, "GRanges")){
        counts = count_dbs_contexts_gr(grl, categories)
        colnames(counts) = "My_sample"
    } else{
        not_gr_or_grl(grl)
    }
    counts = cbind(categories, counts)
    counts[is.na(counts)] = 0
    counts = counts %>% 
        tidyr::unite("muttype_total", REF, ALT) %>% 
        tibble::column_to_rownames("muttype_total") %>% 
        as.matrix()
    
    #counts = tibble::as_tibble(counts)
    #counts$REF = factor(counts$REF, levels = BiocGenerics::unique(counts$REF))
    
    #bases = c("A", "C", "G", "T")
    #bases1 = bases
    #bases_combi = tidyr::crossing(bases, bases1)
    #counts$ALT = factor(counts$ALT, levels = stringr::str_c(bases_combi$bases, bases_combi$bases1))
    return(counts)
}



#' Count dbs contexts from a single GRanges object.
#' 
#' @details
#' Counts the number of DBS per COSMIC context from a GRanges object containing DBS mutations.
#' The function is called by count_dbs_contexts
#' 
#' @param gr GRanges object containing DBS mutations in which the context was added with set_dbs_context.
#' @param categories A tibble containing all possible dbs context categories
#' 
#' @return A single column tibble containing the number of dbs per COSMIC context.
#' 
#' @importFrom magrittr %>%
#' 
#' @seealso \code{\link{count_dbs_contexts}}, \code{\link{set_dbs_context}}
#' @family DBS
count_dbs_contexts_gr = function(gr, categories){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    REF = ALT = NULL
    
    context = cbind("REF" = as.vector(gr$REF), "ALT" = as.vector(unlist(gr$ALT)))
    counts = context %>% 
        tibble::as_tibble() %>% 
        dplyr::group_by(REF, ALT) %>% 
        dplyr::summarise(count = dplyr::n())
    counts_full = dplyr::left_join(categories, counts, by = c("REF", "ALT")) %>% 
        dplyr::select(-REF, -ALT)
    return(counts_full)
}