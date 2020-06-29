#' Plot the similarity between a mutation matrix and its reconstructed profile
#'
#' When a reconstructed profile has a cosine similarity of more than 0.95 with
#' original, the reconstructed profile is considered very good.
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param reconstructed A reconstructed mutation count matrix
#' @param y_intercept The y intercept of the plotted horizontal line. Default: 0.95.
#' @param ylims The limits of the y axis. Default: c(0.6, 1)
#'
#' @return A ggplot figure
#' @export
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Create figure
#' plot_original_vs_reconstructed(mut_mat, nmf_res$reconstructed)
#'
#' ## You can also use the results of signature refitting.
#' ## Here we load some data as an example
#' fit_res <- readRDS(system.file("states/snv_refit.rds",
#'   package = "MutationalPatterns"
#' ))
#' plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed)
#'
#' ## You can also change the height of the horizontal line
#' plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, y_intercept = 0.90)
#'
#' ## It's also possible to change the limits of the y axis
#' plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, ylims = c(0, 1))
plot_original_vs_reconstructed <- function(mut_matrix, reconstructed, y_intercept = 0.95, ylims = c(0.6, 1)) {

  # Determine cosine similarities
  cos_sim_all <- cos_sim_matrix(mut_matrix, reconstructed)
  cos_sim <- as.data.frame(diag(cos_sim_all))

  # Set names
  colnames(cos_sim) <- "cos_sim"
  cos_sim$sample <- row.names(cos_sim_all)

  # Create plot
  plot <- ggplot(cos_sim, aes(y = cos_sim, x = factor(sample, levels = sample))) +
    geom_bar(stat = "identity", fill = "skyblue3") +
    geom_hline(aes(yintercept = y_intercept)) +
    coord_cartesian(ylim = ylims, expand = FALSE) +
    labs(y = "Cosine similarity\n original VS reconstructed", x = "") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  return(plot)
}
