#' Signature clustering function
#'
#' Hierarchical clustering of signatures based on cosine similarity
#'
#' @param signatures Matrix with 96 trinucleotides (rows) and any number of
#' signatures (columns)
#' @param method     The agglomeration method to be used for hierarchical
#' clustering. This should be one of "ward.D", "ward.D2", "single", "complete",
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC). Default = "complete".
#' @return hclust object
#'
#' @examples
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#'
#' ## Hierarchically cluster the cancer signatures based on cosine similarity
#' hclust_signatures <- cluster_signatures(signatures)
#'
#' ## Plot dendrogram
#' plot(hclust_signatures)
#' @seealso
#' \code{\link{plot_contribution_heatmap}}
#'
#' @export

cluster_signatures <- function(signatures, method = "complete") {
  # construct cosine similarity matrix
  sim <- cos_sim_matrix(signatures, signatures)
  # transform to distance
  dist <- as.dist(1 - sim)
  # perform hierarchical clustering
  hc_sig_cos <- hclust(dist, method = method)
  return(hc_sig_cos)
}
