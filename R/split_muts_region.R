#' Split GRangesList or GRanges based on a list of regions.
#'
#' A GRangesList or GRanges object containing variants is split based on a list of regions.
#' This list can be either a GRangesList or a GRanges object.
#' The result is a GRangesList where each element contains the variants of one sample from one region.
#' Variant that are not in any of the provided region are put in a list of 'other'.
#'
#' @param vcf_list GRangesList or GRanges object
#' @param ranges_grl GRangesList or GRanges object containing regions of interest
#' @param include_other Boolean. Whether or not to include a "Other" region
#' containing mutations that aren't in any other region.
#'
#' @return GRangesList
#' @export
#' @family genomic_regions
#' @importFrom magrittr %>%
#' @examples
#'
#' ## Read in some existing genomic regions.
#' ## See the 'genomic_distribution()' example for how we obtained the
#' ## following data:
#' CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#' promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#' flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Combine the regions into a single GRangesList
#' regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
#'
#' names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
#'
#' ## Read in some variants.
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Split muts based on the supplied regions
#' split_muts_region(grl, regions)
#'
#' ## Don't include muts outside of the supplied regions
#' split_muts_region(grl, regions, include_other = FALSE)
split_muts_region <- function(vcf_list, ranges_grl, include_other = TRUE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  # Check arguments
  if (inherits(vcf_list, "CompressedGRangesList")) {
    if (is.null(names(vcf_list))) {
      stop("Please set sample names (without dots) for the vcf_list with `names(vcf_list) = my_names`", call. = FALSE)
    }
    if (any(stringr::str_detect(names(vcf_list), "\\."))) {
      stop("The sample names of the vcf_list should not contain dots. Please fix them with `names(vcf_list) = my_names`", call. = FALSE)
    }
  }
  if (any(stringr::str_detect(names(ranges_grl), "\\."))) {
    stop("The names of the ranges_grl should not contain dots. Please fix them with `names(ranges_grl) = my_names`", call. = FALSE)
  }


  # Turn grl into list.
  if (inherits(vcf_list, "CompressedGRangesList")) {
    vcf_list <- as.list(vcf_list)
  }

  # Get muttype per sample
  if (inherits(vcf_list, "list")) {
    grl <- purrr::map(vcf_list, function(x) .split_muts_region_gr(x, ranges_grl, include_other)) %>%
      purrr::map(as.list) %>% # Create a list of lists. Outer layer: samples. Inner layer: regions.
      do.call(c, .) %>%
      GenomicRanges::GRangesList()
    return(grl)
  } else if (inherits(vcf_list, "GRanges")) {
    grl <- .split_muts_region_gr(vcf_list, ranges_grl, include_other)
    return(grl)
  } else {
    .not_gr_or_grl(vcf_list)
  }
}

#' Split a GRanges object based on a list of regions.
#'
#' A GRanges object containing variants is split based on a list of regions.
#' This list can be either a GRangesList or a GRanges object.
#' The result is a GRangesList where each element contains the variants of one sample from one region.
#' Variant that are not in any of the provided region are put in a list of 'other'.
#'
#' @param gr GRanges object containing variants
#' @param ranges_grl GRangesList or GRanges object containing regions of interest
#' @param include_other Boolean. Whether or not to include a "Other" region
#' containing mutations that aren't in any other region.
#' @noRd
#' @return GRangesList
#'
.split_muts_region_gr <- function(gr, ranges_grl, include_other) {

  # Get the muts overlapping the ranges_grl
  if (inherits(ranges_grl, "CompressedGRangesList")) {
    ranges_list <- as.list(ranges_grl)
    gr_sub_l <- purrr::map(ranges_list, ~ .get_muts_region(gr, .x))
  } else if (inherits(ranges_grl, "GRanges")) {
    gr_sub <- .get_muts_region(gr, ranges_grl)
    gr_sub_l <- list("In_region" = gr_sub)
  } else {
    .not_gr_or_grl(ranges_grl)
  }


  if (include_other) {
    # Get the other muts
    hits <- GenomicRanges::findOverlaps(gr, ranges_grl)
    gr_other <- gr[-S4Vectors::queryHits(hits)]
    gr_sub_l[["Other"]] <- gr_other
  }

  # Combine data
  grl <- GenomicRanges::GRangesList(gr_sub_l)
  return(grl)
}

#' Split a GRanges object based on a list of regions.
#'
#' A GRanges object containing variants is split based on a list of regions.
#' This list can be either a GRangesList or a GRanges object.
#' The result is a GRangesList where each element contains the variants of one sample from one region.
#' Variant that are not in any of the provided region are put in a list of 'other'.
#'
#' @param gr GRanges object containing variants
#' @param ranges GRanges object containing regions of interest
#' @noRd
#' @return GRanges object
#'
.get_muts_region <- function(gr, ranges) {
  hits <- GenomicRanges::findOverlaps(gr, ranges)
  gr_sub <- gr[S4Vectors::queryHits(hits)]
  return(gr_sub)
}
