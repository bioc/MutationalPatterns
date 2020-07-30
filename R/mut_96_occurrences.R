#' Count 96 trinucleotide mutation occurrences
#'
#'  @details
#'  This function is called by mut_matrix. It calculates the 96 trinucleotide context for all variants
#'  and then splits these per GRanges (samples). It then calculates how often each 96 trinucleotide context occurs.
#'
#'
#' @param type_context result from type_context function
#' @param gr_sizes A vector indicating the number of variants per GRanges
#' @return Mutation matrix with 96 trinucleotide mutation occurrences
#'
#' @importFrom magrittr %>%
mut_96_occurrences <- function(type_context, gr_sizes) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  categories <- count <- NULL

  # Determine nr of bases
  nr_bases <- nchar(type_context$context[[1]])
  middle_base <- ceiling(nr_bases / 2)


  # Determine all possible contexts
  bases_left <- c("A", "C", "G", "T")
  bases_right <- c("A", "C", "G", "T")
  base_subs <- c("[C>A]", "[C>G]", "[C>T]", "[T>A]", "[T>C]", "[T>G]")

  # Loop over each base substitution
  full_context_poss <- vector("list", length(base_subs))
  for (i in seq_along(base_subs)) {
    sub <- base_subs[[i]]
    sub_context <- sub
    # Repeatedly add bases left and right
    for (j in seq_len(middle_base - 1)) {
      combi_tb <- tidyr::crossing(bases_left, sub_context, bases_right)
      sub_context <- paste0(combi_tb$bases_left, combi_tb$sub_context, bases_right)
    }
    full_context_poss[[i]] <- sub_context
  }
  full_context_poss <- do.call(c, full_context_poss)


  # Determine 96 context for all variants
  full_context <- stringr::str_c(
    substr(type_context$context, 1, middle_base - 1),
    "[",
    type_context$types,
    "]",
    substr(type_context$context, middle_base + 1, nr_bases)
  ) %>%
    factor(levels = full_context_poss)

  # Set names if they are not yet present
  if (is.null(names(gr_sizes))) {
    names(gr_sizes) <- seq_along(gr_sizes)
  }

  # Create vector describing the sample of each variant
  sample_vector <- rep(names(gr_sizes), gr_sizes) %>%
    factor(levels = names(gr_sizes))

  # Count the mutations per type and per sample
  counts <- tibble::tibble("categories" = full_context, "sample" = sample_vector) %>%
    dplyr::filter(!is.na(categories)) %>%
    dplyr::group_by(categories, sample, .drop = FALSE) %>%
    dplyr::summarise(count = dplyr::n())

  # Transform the data into a mutation matrix
  counts <- tidyr::spread(counts, key = sample, value = count, fill = 0)
  unnecesary_cols <- which(colnames(counts) == "<NA>")
  mut_mat <- as.matrix(counts[, -c(1, unnecesary_cols)])
  rownames(mut_mat) <- counts$categories
  return(mut_mat)
}
