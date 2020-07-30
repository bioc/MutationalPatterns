#' Count 192 trinucleotide mutation occurrences
#'
#'  @details
#'  This function is called by mut_matrix_stranded.
#'  The 192 trinucleotide context is the 96 trinucleotide context combined with the strands.
#'  This function calculates the 192 trinucleotide context for all variants.
#'  and then splits these per GRanges (samples). It then calculates how often each 192 trinucleotide context occurs.
#'
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand,
#' and "-" for unknown
#' @param gr_sizes A vector indicating the number of variants per GRanges
#'
#' @importFrom magrittr %>%
#'
#' @return Mutation matrix with 192 mutation occurrences and 96 trinucleotides
#' for two strands

mut_192_occurrences <- function(type_context, strand, gr_sizes) {
  # get possible strand values
  values <- levels(strand)

  idx1 <- which(strand == values[1])
  idx2 <- which(strand == values[2])

  # get type context for both vcf subsets
  type_context_1 <- purrr::map(type_context, function(x) x[idx1])
  type_context_2 <- purrr::map(type_context, function(x) x[idx2])

  # Subset the gr_sizes.
  sample_vector <- rep(names(gr_sizes), gr_sizes) %>%
    factor(levels = names(gr_sizes))
  table_vector_1 <- sample_vector[idx1] %>%
    table()
  gr_sizes_1 <- as.vector(table_vector_1)
  names(gr_sizes_1) <- names(table_vector_1)
  table_vector_2 <- sample_vector[idx2] %>%
    table()
  gr_sizes_2 <- as.vector(table_vector_2)
  names(gr_sizes_2) <- names(table_vector_2)

  # make 96-trinucleotide count vector per set
  mut_mat_1 <- mut_96_occurrences(type_context_1, gr_sizes_1)
  mut_mat_2 <- mut_96_occurrences(type_context_2, gr_sizes_2)

  # add names
  names_1 <- paste(rownames(mut_mat_1), values[1], sep = "-")
  names_2 <- paste(rownames(mut_mat_2), values[2], sep = "-")

  # combine matrixes
  mut_mat <- rbind(mut_mat_1, mut_mat_2)
  rownames(mut_mat) <- c(names_1, names_2)

  # Reorder for backwards compatibility
  reorder_i <- purrr::map2(
    seq(1, nrow(mut_mat) / 2),
    seq(
      nrow(mut_mat) / 2 + 1,
      nrow(mut_mat)
    ),
    c
  ) %>%
    unlist()
  mut_mat <- mut_mat[reorder_i, , drop = FALSE]

  return(mut_mat)
}
