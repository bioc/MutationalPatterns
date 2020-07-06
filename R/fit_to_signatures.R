#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#'
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving the nonnegative least-squares
#' constraints problem.
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param signatures Signature matrix (dimensions: x mutation types
#' X n signatures)
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @importFrom pracma lsqnonneg
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get signatures
#' signatures <- get_known_signatures()
#'
#' ## Perform the fitting
#' fit_res <- fit_to_signatures(mut_mat, signatures)
#'
#' ## This will also work for indels and dbs.
#' ## An example is given for indels
#'
#' ## Get The indel counts
#' ## See 'count_indel_contexts()' for more info on how to do this.
#' indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get signatures
#' signatures <- get_known_signatures("indel")
#'
#' fit_to_signatures(indel_counts, signatures)
#' @seealso \code{\link{mut_matrix}},\code{\link{fit_to_signatures_strict}},\code{\link{fit_to_signatures_bootstrapped}}
#'
#' @export

fit_to_signatures <- function(mut_matrix, signatures) {
  # make sure dimensions of input matrix are correct
  if (dim(mut_matrix)[1] != dim(signatures)[1]) {
    stop(paste(
      "Mutation matrix and signatures input have",
      "different number of mutational features"
    ))
  }

  n_features <- dim(mut_matrix)[1]
  n_samples <- dim(mut_matrix)[2]
  n_signatures <- dim(signatures)[2]
  lsq_contribution <- matrix(NA, nrow = n_signatures, ncol = n_samples)
  lsq_reconstructed <- matrix(NA, nrow = n_features, ncol = n_samples)

  # Process each sample
  for (i in seq_len(ncol(mut_matrix)))
  {
    y <- mut_matrix[, i]
    lsq <- lsqnonneg(signatures, y)
    lsq_contribution[, i] <- lsq$x
    lsq_reconstructed[, i] <- signatures %*% as.matrix(lsq$x)
  }

  # Add row and col names
  sample_names <- colnames(mut_matrix)
  signature_names <- colnames(signatures)
  mut_type_names <- rownames(signatures)

  colnames(lsq_contribution) <- sample_names
  rownames(lsq_contribution) <- signature_names

  colnames(lsq_reconstructed) <- sample_names
  rownames(lsq_reconstructed) <- mut_type_names

  res <- list(lsq_contribution, lsq_reconstructed)
  names(res) <- c("contribution", "reconstructed")

  return(res)
}
