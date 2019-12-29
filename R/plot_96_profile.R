#' Plot 96 trinucleotide profile
#'
#' Plot relative contribution of 96 trinucleotides      
#' @param mut_matrix 96 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors 6 value color vector
#' @param condensed More condensed plotting format. Default = F.
#' @return 96 trinucleotide profile plot
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## mutation matrix with transcriptional strand information:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Plot the 96-profile of three samples
#' plot_96_profile(mut_mat[,c(1,4,7)])
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

plot_96_profile = function(mut_matrix, colors, ymax = 0.2, condensed = FALSE)
{
    warning("Function will be deprecated. Use 'plot_profiles' instead")
    plot = plot_profiles(mut_matrix = mut_matrix, type = "snv", 
                         colors = colors, ymax = c("snv"=ymax), 
                         condensed = condensed)
    return(plot)
}
