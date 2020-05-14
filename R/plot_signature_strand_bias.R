#' Plot signature strand bias
#' 
#' Plot strand bias per mutation type for each signature.
#' 
#' @param signatures_strand_bias Signature matrix with 192 features
#' @return Barplot
#'
#' @import ggplot2
#' @importFrom magrittr %>% 
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
#'
#' nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
#'                                         package="MutationalPatterns"))
#'
#' ## Provide column names for the plot.
#' colnames(nmf_res_strand$signatures) = c("Signature A", "Signature B")
#'
#' plot_signature_strand_bias(nmf_res_strand$signatures) 
#'
#' @seealso
#' \code{link{extract_signatures}},
#' \code{link{mut_matrix}}
#'
#' @export

plot_signature_strand_bias = function(signatures_strand_bias)
{
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    Signature = type = ratio = transcribed = untranscribed = NULL
    `Group.1` = `Group.2` = significant = pval = NULL

    # check if there are 192 features in the signatures
    if (dim(signatures_strand_bias)[1] != 192)
        stop(paste("Input signature matrix does not have 192 features (96",
                    "trinucleotide * 2 strands)."))

    # aggregate by strand and type
    sum_per_type = stats::aggregate(signatures_strand_bias,
                                    by=list(STRAND, SUBSTITUTIONS_192),
                                    FUN=sum)

    #Make tibble longer
    stats_per_type = sum_per_type %>% 
        dplyr::rename(strand = `Group.1`, type = `Group.2`) %>% 
        dplyr::mutate(strand = ifelse(strand == "T", "transcribed", "untranscribed")) %>% 
        tidyr::pivot_longer(c(-strand, -type), names_to = "Signature") %>%  #Combine signature columns
        tidyr::pivot_wider(names_from = strand) %>% #Split transcribed/untranscribed column for binom test.
        dplyr::mutate(observed = as.integer(transcribed),
                      size = as.integer(transcribed + untranscribed),
                      ratio = transcribed / untranscribed)
    
    # Perform binomial test
    stats_per_type$pval = purrr::map_dbl(seq_len(nrow(stats_per_type)), function(i){
        row = stats_per_type[i,]
        binomial_test(0.5, row$size, row$observed)$pval
    })
    #Add stars when significant.
    stats_per_type = dplyr::mutate(stats_per_type, significant = ifelse(pval < 0.05, "*", " "))
   
    # Find maximum y value for plotting
    max = round(max(abs(log2(stats_per_type$ratio))))

    # Plot
    plot = ggplot(stats_per_type,
                    aes(x=type, y=log2(ratio), fill=type)) +
        geom_bar(stat="identity", position="dodge", color="black") +
        scale_y_continuous(limits=c(-max,max)) +
        scale_fill_manual(values=COLORS6) +
        facet_grid(Signature ~ .) +
        ylab("log2(transcribed/untranscribed)") +
        theme_bw() + 
        scale_x_discrete(breaks=NULL) +
        xlab("") +
        geom_text(
            aes(x = type,
                y = log2(ratio),
                label = significant,
                vjust = ifelse(sign(log2(ratio)) > 0, 0.5, 1)), 
            size = 8, position = ggplot2::position_dodge(width=1))

    return(plot)
}
