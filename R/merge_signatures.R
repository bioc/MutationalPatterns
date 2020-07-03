#' Merge signatures based on cosine similarity
#'
#' This function merges signatures based on their cosine similarity.
#' It iteratively merges the two signatures with the highest cosine similarity.
#' Merging is stopped when the maximum cosine similarity is lower than the limit.
#'
#' @param signatures Signature matrix (dimensions: x mutation types
#' X n signatures)
#' @param cos_sim_cutoff Cutoff for cosine similarity. Signatures are merged when their
#' cosine similarity is higher than the limit. Default: 0.8
#' @param merge_char Character used to merge the signature names. This character shouldn't
#' be in the signature names beforehand. Default: ";"
#' @param verbose Verbosity. If TRUE it shows which signatures got merged. Default: TRUE
#'
#' @return Signature matrix (dimensions: x mutation types
#' X n signatures)
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## Merge signatures
#' merge_signatures(signatures)
#'
#'
#' ## Merge signatures using a stricter cutoff
#' merge_signatures(signatures, cos_sim_cutoff = 0.9)
#'
#' ## Merge signatures using a different merging character
#' merge_signatures(signatures, merge_char = "_")
#'
#' ## Merge signatures silently
#' merge_signatures(signatures, verbose = FALSE)
merge_signatures <- function(signatures, cos_sim_cutoff = 0.8, merge_char = ";", verbose = TRUE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  # Validate signature names
  nr_lowercase <- colnames(signatures) %>%
    grepl(merge_char, .) %>%
    sum()
  if (nr_lowercase > 0) {
    stop(paste0("Please remove all ", merge_char, " characters from your signature names."), call. = FALSE)
  }

  # Determine max similarity between signatures
  sim_m <- cos_sim_matrix(signatures, signatures)
  diag(sim_m) <- 0
  max <- max(sim_m)

  # Merge signatures while max similarity is higher than cutoff.
  while (max > cos_sim_cutoff) {

    # Find signatures that need to be merged
    max_index <- order(sim_m, decreasing = TRUE)[1]
    max_loc <- arrayInd(max_index, dim(sim_m), useNames = TRUE)
    sigs_left <- signatures[, -max_loc, drop = FALSE]
    sigs_to_combi <- signatures[, max_loc, drop = FALSE]

    # Signatures that have already been merged and thus exist
    # of multiple signatures are weighted accordingly.
    weights <- sigs_to_combi %>%
      colnames() %>%
      stringr::str_count(merge_char) %>%
      magrittr::add(1)
    combi_sig <- sigs_to_combi %*% diag(weights)

    # Merge signatures
    combi_sig <- combi_sig %>%
      rowSums() %>%
      matrix()
    combi_sig <- combi_sig / sum(weights)
    colnames(combi_sig) <- paste(colnames(sigs_to_combi), collapse = merge_char)

    # Add merged signature to the rest.
    signatures <- cbind(sigs_left, combi_sig)

    # Print which signatures have been merged
    if (verbose) {
      merged_sig_names <- paste0(colnames(sigs_to_combi), collapse = ", ")
      message(paste0("Combined the following two signatures: ", merged_sig_names))
    }

    # Determine max similarity between signatures for next loop
    sim_m <- cos_sim_matrix(signatures, signatures)
    diag(sim_m) <- 0
    max <- max(sim_m)
  }
  return(signatures)
}
