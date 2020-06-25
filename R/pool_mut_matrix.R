#' Pool multiple samples from a mutation matrix together
#'
#' The mutation counts of columns (samples) are added up according to the grouping variable.
#'
#' @param mut_matrix Mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param grouping Grouping variable
#'
#' @return Mutation count matrix (dimensions: x mutation types
#' X n groups)
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#' grouping <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#' pool_mut_mat(mut_mat, grouping)
pool_mut_mat <- function(mut_matrix, grouping) {
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  grouping <- factor(grouping)
  mut_mat_group <- mut_matrix %>%
    t(.) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(factor = grouping) %>%
    dplyr::group_by(factor) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::select(-factor) %>%
    t(.)
  colnames(mut_mat_group) <- levels(grouping)
  return(mut_mat_group)
}
