#' Significance test for strand asymmetry
#'
#' This function performs a two sided Poisson test for the ratio between mutations on 
#' each strand. Multiple testing correction is also performed.
#' 
#' @param strand_occurrences Dataframe with mutation count per strand, result
#' from strand_occurrences()
#' @return Dataframe with poisson test P value for the ratio between the
#' two strands per group per base substitution type.
#' @importFrom magrittr  %>% 
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' ## Perform the strand bias test.
#' strand_counts = strand_occurrences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurrences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_bias_test = function(strand_occurrences){
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    group = type = strand = variable = relative_contribution = no_mutations = p_poisson = NULL

    #Make data long
    df_strand = strand_occurrences %>% 
        dplyr::select(-relative_contribution) %>% 
        tidyr::pivot_wider(names_from = strand, values_from = no_mutations)
        
    #Calculate total and ratio
    df_strand[,"total"] = df_strand[,3] + df_strand[,4]
    df_strand[,"ratio"] = df_strand[,3] / df_strand[,4]
    
    #poisson test for strand ratio
    #Is the same as binom in this scenario. 
    #Since the output uses poisson in the name, we will keep using this.
    df_strand$p_poisson = apply(df_strand, 1, function(x){
        stats::poisson.test(c(as.numeric(x[3]), as.numeric(x[4])), 
                            r=1, alternative = "two.sided")$p.value
        })
    
    #Add significance stars and do multiple testing correction.
    df_strand = df_strand %>% 
        dplyr::mutate(significant = ifelse(p_poisson < 0.05, "*", " "),
                      fdr = stats::p.adjust(p_poisson, method = "fdr"),
                      significant_fdr = ifelse(fdr < 0.1, "*", " "))
    
    return(df_strand)
}
