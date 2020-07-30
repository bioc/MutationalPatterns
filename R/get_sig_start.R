#' Determine the number of significance stars
#'
#' The number of significance stars is determined based on the statistical value
#' and the significance cutoffs.
#'
#' @param val Statistical value. Either a p value or fdr.
#' @param cutoffs Significance cutoffs for the statistical value.
#'
#' @return A vector of significance stars and empty strings (not significant).
#' @noRd
#'
.get_sig_star <- function(val, cutoffs) {

  # Get name of cutoffs argument
  cutoffs_name <- deparse(substitute(cutoffs))

  # Validate cutoff argument
  if (length(cutoffs) > 3) {
    stop(paste0("The length of the ", cutoffs_name, " argument can't be higher than 3."),
      call. = FALSE
    )
  }

  if (!all.equal(cutoffs, sort(cutoffs, decreasing = TRUE))) {
    stop(paste0("The ", cutoffs_name, " argument should be in decreasing order."),
      call. = FALSE
    )
  }

  # Add -Infs to cutoffs if the length is lower than 3.
  # Since a val cant be lower than -Inf, these cutoffs will never be reached.
  cutoffs <- c(cutoffs, rep(-Inf, 3 - length(cutoffs)))


  # Determine significance level
  stars <- dplyr::case_when(
    val < cutoffs[3] ~ "***",
    val < cutoffs[2] ~ "**",
    val < cutoffs[1] ~ "*",
    TRUE ~ ""
  )
  return(stars)
}
