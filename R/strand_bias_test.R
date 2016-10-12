#' Significance test for transcriptional strand asymmetry
#'
#' This function performs a Poisson test for the ratio between mutations on the
#' transcribed and untranscribed strand
#' @param strand_occurences Dataframe with mutation count per strand, result
#' from strand_occurences()
#' @return Dataframe with poisson test P value for the ratio between the
#' transcribed and untrascribed strand per group per base substitution type.
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom plyr .
#' @importFrom plyr ddply
#' @importFrom plyr summarise
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
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurences}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_bias_test = function(strand_occurences)
{
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    variable = NULL
    transcribed = NULL
    untranscribed = NULL

    # statistical test for strand ratio
    # poisson test
    df_strand = reshape2::dcast(melt(strand_occurences),
                                type + group ~ strand,
                                sum,
                                subset = plyr::.(variable == "no_mutations"))

    ## Prevent using 'T' as a variable name to avoid confusion with TRUE.
    ## Rename the columns 'T' and 'U' to 'Ts' and 'Us'.
    colnames(df_strand) <- c("type", "group", "transcribed", "untranscribed")

    df_strand = plyr::ddply(df_strand,
                            c("group", "type", "transcribed", "untranscribed"),
                            plyr::summarise,
                            total = transcribed+untranscribed,
                            ratio = transcribed/untranscribed,
                            p_poisson = poisson.test(c(untranscribed,
                                                        transcribed),
                                                    r=1)$p.value)

    df_strand$significant[df_strand$p_poisson < 0.05] = "*"
    df_strand$significant[df_strand$p_poisson >= 0.05] = " "

    return(df_strand)
}
