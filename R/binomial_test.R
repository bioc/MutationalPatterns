#' Binomial test for enrichment or depletion testing
#'
#' This function performs lower-tail binomial test for depletion and
#' upper-tail test for enrichment
#'
#' @param p Probability of success
#' @param n Number of trials
#' @param x Observed number of successes
#' @param p_cutoffs Significance cutoff for the p value. Default: 0.05
#' @return A data.frame with direction of effect (enrichment/depletion),
#' P-value and significance asterisks
#'
#' @examples
#' binomial_test(0.5, 1200, 543)
#' binomial_test(0.2, 800, 150)
#' @export

binomial_test <- function(p, n, x, p_cutoffs = 0.05) {
  # Calculate expected number of successes
  expected <- p * n

  # Handle depletion
  if (x < expected) {
    # do lower tail test
    pval <- stats::pbinom(x, n, p, lower.tail = TRUE)
    effect <- "depletion"
  }

  # Handle enrichment
  else {
    # do upper tail test
    pval <- stats::pbinom(x - 1, n, p, lower.tail = FALSE)
    effect <- "enrichment"
  }

  # make test two sided.
  pval <- 2 * min(pval, 1 - pval)

  # Add significance asteriks
  significant <- .get_sig_star(pval, p_cutoffs)

  res <- data.frame("effect" = factor(effect), pval, significant)
  return(res)
}
