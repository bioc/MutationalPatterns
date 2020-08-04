#' Count occurrences per base substitution type and strand
#'
#' For each base substitution type and strand the total number
#' of mutations and the relative contribution within a group is returned.
#'
#' @param mut_mat_s 192 feature mutation count matrix, result from
#' 'mut_matrix_stranded()'
#' @param by Character vector with grouping info, optional
#'
#' @return A data.frame with the total number of mutations and relative
#' contribution within group per base substitution type and strand
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c(
#'   "colon", "colon", "colon",
#'   "intestine", "intestine", "intestine",
#'   "liver", "liver", "liver"
#' )
#'
#' strand_counts <- strand_occurrences(mut_mat_s, by = tissue)
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{plot_strand}},
#' \code{\link{plot_strand_bias}}
#'
#' @export

strand_occurrences <- function(mut_mat_s, by = NA) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  type_strand <- no_mutations <- group <- NULL

  # Make data long
  tb_per_sample <- mut_mat_s %>%
    as.data.frame() %>%
    tibble::rownames_to_column("type_strand") %>%
    tidyr::pivot_longer(-type_strand, values_to = "no_mutations", names_to = "sample") %>%
    tidyr::separate(type_strand, into = c("type", "strand"), sep = "-") %>%
    dplyr::mutate(type = stringr::str_replace(type, ".*\\[(.*)\\].*", "\\1"))

  # If grouping variable not provided, set to "all"
  if (.is_na(by)) {
    by <- "all"
  }
  
  # Add grouping info
  tb_by <- tibble::tibble(
    "sample" = unique(tb_per_sample$sample),
    "by" = by
  )
  tb_per_sample <- tb_per_sample %>%
    dplyr::left_join(tb_by, by = "sample")

  # Summarise per group
  tb <- tb_per_sample %>%
    dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, type, strand) %>% # Summarise the number of mutations per group
    dplyr::summarise(no_mutations = sum(no_mutations)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(by) %>% # Calculate relative per group. NOT per sample.
    dplyr::mutate(relative_contribution = no_mutations / sum(no_mutations)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(group = by)

  return(tb)
}
