#' Plot strand bias per base substitution type
#'
#' For each base substitution type and transcriptional strand the total
#' number of mutations and the relative contribution within a group is
#' returned.
#'
#' @param strand_bias data.frame, result from strand_bias function
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#'
#' @import ggplot2
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' ## Perform the strand bias test.
#' strand_counts = strand_occurences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' ## Plot the strand bias.
#' plot_strand_bias(strand_bias)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurences}},
#' \code{\link{strand_bias_test}}
#' \code{\link{plot_strand}}
#'
#' @export

plot_strand_bias = function(strand_bias, colors)
{
    # if colors parameter not provided, set to default colors
    if (missing(colors))
        colors=COLORS6

    # determine max y value for plotting
    # = log2 ratio with pseudo counts of 0.1
    log2_ratio = log2(  (strand_bias$transcribed+0.1) /
                        (strand_bias$untranscribed+0.1))

    # max yvalue for plotting plus
    max = round(max(abs(log2_ratio)), digits = 1) + 0.1

    type = NULL
    ratio = NULL
    significant = NULL
    transcribed = NULL
    untranscribed = NULL

    # add label for infinite values
    label2 = log2(strand_bias$ratio)
    select = which(is.finite(label2))
    label2[select] = " "

    # plot strand bias with poisson test results
    plot = ggplot(strand_bias, aes( x = type,
                                    y = log2((transcribed+0.1) /
                                                (untranscribed+0.1)),
                                    fill = type)) +
        scale_fill_manual(values = COLORS6) +
        geom_bar(colour = "black", stat ="identity", position = "identity") +
        scale_y_continuous(limits = c(-max, max)) +
        geom_text(
            aes(x = type,
                y = log2((transcribed+0.1) / (untranscribed+0.1)),
                ymax = log2((transcribed+0.1) / (untranscribed+0.1)), 
                label = significant,
                vjust = ifelse(sign(log2((transcribed+0.1) /
                                            (untranscribed+0.1))) > 0, 0.5, 1)),
            size = 8,
            position = ggplot2::position_dodge(width = 1)) +
        # geom_text(
        #     aes(x = type,
        #         y = log2((transcribed+0.1) / (untranscribed+0.1)),
        #         ymax = log2((transcribed+0.1) / (untranscribed+0.1)), 
        #         label = label2,
        #         vjust = ifelse(sign(log2((transcribed+0.1) /
        #                                  (untranscribed+0.1))) > 0, 0.5, 1)),
        #     size = 3,
        #     position = ggplot2::position_dodge(width = 1)) +
        facet_grid(. ~ group) +
        theme_bw()  +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                legend.title = element_blank()) +
        xlab("") + 
        ylab("log2(transcribed/untranscribed)") +
        scale_x_discrete(breaks=NULL)

    return(plot)
}
