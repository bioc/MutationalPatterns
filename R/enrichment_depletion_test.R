#' Test for enrichment or depletion of mutations in genomic regions
#'
#' This function aggregates mutations per group (optional) and performs an
#' enrichment depletion test.
#'
#' @param x data.frame result from genomic_distribution()
#' @param by Optional grouping variable, e.g. tissue type
#' @param p_cutoffs Significance cutoff for the p value. Default: 0.05
#' @param fdr_cutoffs Significance cutoff for the fdr. Default: 0.1
#' @return data.frame with the observed and expected number of mutations per
#' genomic region per group (by) or sample
#'
#'
#' @examples
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' distr <- readRDS(system.file("states/distr_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
#'
#' ## Perform the enrichment/depletion test by tissue type.
#' distr_test <- enrichment_depletion_test(distr, by = tissue)
#'
#' ## Or without specifying the 'by' parameter, to pool all samples.
#' distr_single_sample <- enrichment_depletion_test(distr)
#'
#' ## Use different significance cutoffs for the pvalue and fdr
#' distr_strict <- enrichment_depletion_test(distr,
#'   by = tissue,
#'   p_cutoffs = 0.01, fdr_cutoffs = 0.05
#' )
#'
#' ## Use multiple (max 3) significance cutoffs.
#' ## This will vary the number of significance stars.
#' distr_multistars <- enrichment_depletion_test(distr,
#'   by = tissue,
#'   p_cutoffs = c(0.05, 0.01, 0.005),
#'   fdr_cutoffs = c(0.1, 0.05, 0.01)
#' )
#' @seealso
#' \code{\link{genomic_distribution}},
#' \code{\link{plot_enrichment_depletion}}
#'
#' @export

enrichment_depletion_test <- function(x, by = NA,
                                      p_cutoffs = 0.05,
                                      fdr_cutoffs = 0.1) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  pval <- fdr <- sample <- prob <- expected <- region <- n_muts <- NULL
  surveyed_length <- surveyed_region_length <- observed <- NULL
  
  # If grouping variable not provided, set to "all"
  if (.is_na(by)) {
    by <- "all"
  }
  
  # Add grouping info
  tb_by <- tibble::tibble(
    "sample" = unique(x$sample),
    "by" = by
  )
  tb_per_sample <- x %>%
    dplyr::left_join(tb_by, by = "sample")
  
  # Summarize based on the grouping.
  tb <- tb_per_sample %>% 
    dplyr::select(-prob, -expected) %>% 
    dplyr::mutate(region = factor(region, levels = unique(region)),
                                  by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(region, by) %>% 
    dplyr::summarise(n_muts = sum(n_muts), 
                     surveyed_length = sum(surveyed_length), 
                     surveyed_region_length = sum(surveyed_region_length),
                     observed = sum(observed)) %>% 
    dplyr::ungroup()
  
  
  # Calculate probability and expected number of mutations
  tb$prob <- tb$n_muts / tb$surveyed_length
  tb$expected <- tb$prob * tb$surveyed_region_length

  # Perform enrichment/depletion test for each row
  nr_muts <- nrow(tb)
  tb2 <- vector("list", nr_muts)
  for (i in seq_len(nr_muts)) {
    x <- tb[i, ]
    tb2[[i]] <- binomial_test(
      x$prob,
      x$surveyed_region_length,
      x$observed,
      p_cutoffs
    )
  }
  tb2 <- do.call(rbind, tb2)

  # Combine results into one data frame
  df <- cbind(tb, tb2)

  # Calculate fdr
  df <- dplyr::mutate(df,
    fdr = stats::p.adjust(pval, method = "fdr"),
    significant_fdr = .get_sig_star(fdr, fdr_cutoffs)
  )

  return(df)
}
