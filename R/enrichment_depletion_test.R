#' Test for enrichment or depletion of mutations in genomic regions
#'
#' This function aggregates mutations per group (optional) and performs an
#' enrichment depletion test.
#'
#' @param x data.frame result from genomic_distribution() 
#' @param by Optional grouping variable, e.g. tissue type
#' @param p_cutoffs Significance cutoff for the p value. Default: 0.05
#' @param fdr_cutoffs Significance cutoff for the fdr. Default: 0.1
#' @return data.frame with the observed and expected number of mutations per
#' genomic region per group (by) or sample
#'
#'
#' @examples
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' distr <- readRDS(system.file("states/distr_data.rds",
#'                     package="MutationalPatterns"))
#' 
#' tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#'
#' ## Perform the enrichment/depletion test by tissue type.
#' distr_test <- enrichment_depletion_test(distr, by = tissue)
#'
#' ## Or without specifying the 'by' parameter.
#' distr_single_sample <- enrichment_depletion_test(distr)
#' 
#' ## Use different significance cutoffs for the pvalue and fdr
#' distr_strict <- enrichment_depletion_test(distr, by = tissue, 
#'                                          p_cutoffs = 0.01, fdr_cutoffs = 0.05)
#'                                          
#' ## Use multiple (max 3) significance cutoffs.
#' ## This will vary the number of significance stars.
#' distr_multistars <- enrichment_depletion_test(distr, by = tissue, 
#'                                          p_cutoffs = c(0.05, 0.01, 0.005), 
#'                                          fdr_cutoffs = c(0.1, 0.05, 0.01))
#' @seealso
#' \code{\link{genomic_distribution}},
#' \code{\link{plot_enrichment_depletion}}
#'
#' @export

enrichment_depletion_test = function(x, by = c(), 
                                     p_cutoffs = 0.05, 
                                     fdr_cutoffs = 0.1){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    pval = fdr = NULL
    
    # Handle the 'by' parameter when necessary by aggregating x
    if (length(by) > 0){
        x$by = by
        # Sum the columns while aggregating rows based on unique values
        # in 'by' and 'region'.
        res2 = stats::aggregate(cbind(n_muts,
                                        surveyed_length,
                                        surveyed_region_length,
                                        observed) ~ by + region,
                                data = x, sum)
    } else{
        res2 = x
        # In this case, the 'by' variable is 'sample' variable.
        res2$by = res2$sample
        # Select output columns
        res2 = res2[,c(9,1,3,4,6,8)]
    }

    # Calculate probability and expected number of mutations
    res2$prob = res2$n_muts / res2$surveyed_length
    res2$expected = res2$prob * res2$surveyed_region_length

    # Perform enrichment/depletion test for each row
    nr_muts = nrow(res2)
    res3 = vector("list", nr_muts)
    for(i in seq_len(nr_muts)){
        x = res2[i,]
        res3[[i]] = binomial_test(x$prob,
                                  x$surveyed_region_length,
                                  x$observed,
                                  p_cutoffs)
    }
    res3 = do.call(rbind, res3)

    # Combine results into one data frame
    df = cbind(res2, res3)
    
    #Calculate fdr
    df = dplyr::mutate(df, fdr = stats::p.adjust(pval, method = "fdr"),
                  significant_fdr = get_sig_star(fdr, fdr_cutoffs))
    
    return(df)
}
