#' Convert tissue specific signature exposures to reference
#'
#' This function converts tissue specific signature contributions into
#' reference signature contributions. This works on SNV signatures from SIGNAL.
#' It uses a conversion matrix to do the conversion.
#' The output can include possible artifact signatures.
#'
#' @param fit_res Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @return The input fit_res, but with converted signature contributions.
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Get tissue specific signatures
#' signatures <- get_known_signatures(source = "SIGNAL", sig_type = "tissue", tissue_type = "Skin")
#'
#' ## Fit tissue specific signatures
#' fit_res <- fit_to_signatures(mut_mat, signatures)
#'
#' ## Convert the tissue specific signatures exposures to reference
#' fit_res <- convert_sigs_to_ref(fit_res)
convert_sigs_to_ref <- function(fit_res) {

  # Get contribution
  contri <- fit_res$contribution

  # Determine convertion matrix filename
  fname_matrix <- file.path("extdata", "signatures", "SIGNAL_conversion_matrix.txt")
  fname_matrix <- system.file(fname_matrix, package = "MutationalPatterns")

  # Read conversion matrix
  conv_m <- read.table(fname_matrix,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    dec = ",",
    check.names = FALSE
  ) %>%
    tibble::column_to_rownames("Tissue_sig") %>%
    as.matrix()

  # Check that the tissue specific signature names are all in the conversion matrix.
  if (sum(!rownames(contri) %in% rownames(conv_m))) {
    stop(paste0(
      "The signature names of the fit_res don't match that of ",
      "the conversion matrix.\n You have to use tissue specific SNV ",
      "signatures from SIGNAL."
    ), call. = FALSE)
  }

  # Remove signatures, that weren't used, from the conversion matrix.
  conv_m <- conv_m[rownames(conv_m) %in% rownames(contri), , drop = FALSE]


  # Convert signatures to reference.
  fit_res$contribution <- t(conv_m) %*% contri

  return(fit_res)
}
