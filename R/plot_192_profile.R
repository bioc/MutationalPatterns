#' Plot 192 trinucleotide profile
#'
#' Plot relative contribution of 192 trinucleotides      
#' @param mut_matrix 192 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors 6 value color vector
#' @param condensed More condensed plotting format. Default = F.
#' @return 192 trinucleotide profile plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## mutation matrix with transcriptional strand information:
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Extract the signatures.
#' ## This is a computationally intensive task, so we load a precomputed
#' ## version instead.
#' # nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
#' nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Optionally, provide signature names
#' colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")
#'
#' ## Generate the plot
#' plot_192_profile(nmf_res_strand$signatures)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{extract_signatures}}
#'
#' @export

plot_192_profile = function(mut_matrix, colors, ymax = 0.2, condensed = FALSE)
{
    if (grepl("transcribed", rownames(mut_matrix)[1]))
      mode = "transcription"
    else if (grepl("left", rownames(mut_matrix)[1]))
      mode = "replication"

    warning("Function will be deprecated. Use 'plot_strand_profiles' instead")
    plot = plot_strand_profiles(mut_matrix = mut_matrix, type = "snv",
                         colors = colors, ymax = c("snv"=ymax), mode = mode,
                         condensed = condensed)
    return(plot)
}
