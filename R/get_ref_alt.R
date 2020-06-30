#' get REF column from GRanges
#'
#' Retreives the REF column from a GRanges object.
#' This can be spelled as REF, Ref or ref.
#'
#' @param gr GRanges object
#'
#' @return DNAStringSet
#' @noRd
#'
.get_ref <- function(gr) {
  gr_cols <- colnames(S4Vectors::mcols(gr))
  if ("REF" %in% gr_cols) {
    ref <- gr$REF
  } else if ("ref" %in% gr_cols) {
    ref <- gr$ref
  } else if ("Ref" %in% gr_cols) {
    ref <- gr$Ref
  } else {
    stop("Some of your data is missing a REF column.", call. = FALSE)
    ref <- Biostrings::DNAStringSet()
  }
  return(ref)
}

#' get ALT column from GRanges
#'
#' Retreives the ALT column from a GRanges object.
#' This can be spelled as ALT, Alt or alt
#'
#' @param gr GRanges object
#'
#' @return DNAStringSetList
#' @noRd
#'
.get_alt <- function(gr) {
  gr_cols <- colnames(S4Vectors::mcols(gr))
  if ("ALT" %in% gr_cols) {
    alt <- gr$ALT
  } else if ("alt" %in% gr_cols) {
    alt <- gr$alt
  } else if ("Alt" %in% gr_cols) {
    alt <- gr$Alt
  } else {
    stop("Some of your data is missing a ALT column.", call. = FALSE)
    alt <- Biostrings::DNAStringSetList()
  }
  return(alt)
}
