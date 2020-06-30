#' Find overlaps between mutations and a genomic region.
#'
#' Function finds the number of mutations that reside in genomic region and
#' takes surveyed area of genome into account.
#'
#' @param vcf_list GRangesList or GRanges object.
#' @param surveyed_list A GRangesList or a list with GRanges of regions of
#' the genome that have been surveyed (e.g. determined using GATK CallableLoci).
#' @param region_list A GRangesList or a list with GRanges objects containing
#' locations of genomic regions.
#'
#' @return A data.frame containing the number observed and number of expected
#' mutations in each genomic region.
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Use biomaRt to obtain data.
#' ## We can query the BioMart database, but this may take a long time,
#' ## so we take some shortcuts by loading the results from our
#' ## examples. The corresponding code for downloading the data can be
#' ## found above the command we run.
#'
#' # mart="ensemble"
#' # library(biomaRt)
#'
#' # regulatory <- useEnsembl(biomart="regulation",
#' #                          dataset="hsapiens_regulatory_feature",
#' #                          GRCh = 37)
#' regulatory <- readRDS(system.file("states/regulatory_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Download the regulatory CTCF binding sites and convert them to
#' ## a GRanges object.
#' # CTCF <- getBM(attributes = c('chromosome_name',
#' #                             'chromosome_start',
#' #                             'chromosome_end',
#' #                             'feature_type_name'),
#' #              filters = "regulatory_feature_type_name",
#' #              values = "CTCF Binding Site",
#' #              mart = regulatory)
#' #
#' # CTCF_g <- reduce(GRanges(CTCF$chromosome_name,
#' #                 IRanges(CTCF$chromosome_start,
#' #                 CTCF$chromosome_end)))
#'
#' CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Download the promoter regions and conver them to a GRanges object.
#' # promoter = getBM(attributes = c('chromosome_name', 'chromosome_start',
#' #                                 'chromosome_end', 'feature_type_name'),
#' #                  filters = "regulatory_feature_type_name",
#' #                  values = "Promoter",
#' #                  mart = regulatory)
#' # promoter_g = reduce(GRanges(promoter$chromosome_name,
#' #                     IRanges(promoter$chromosome_start,
#' #                             promoter$chromosome_end)))
#'
#' promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' # flanking = getBM(attributes = c('chromosome_name',
#' #                                 'chromosome_start',
#' #                                 'chromosome_end',
#' #                                 'feature_type_name'),
#' #                  filters = "regulatory_feature_type_name",
#' #                  values = "Promoter Flanking Region",
#' #                  mart = regulatory)
#' # flanking_g = reduce(GRanges(
#' #                        flanking$chromosome_name,
#' #                        IRanges(flanking$chromosome_start,
#' #                        flanking$chromosome_end)))
#'
#' flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
#'
#' names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
#'
#' # Use a naming standard consistently.
#' seqlevelsStyle(regions) <- "UCSC"
#'
#' ## Get the filename with surveyed/callable regions
#' surveyed_file <- system.file("extdata/callableloci-sample.bed",
#'   package = "MutationalPatterns"
#' )
#'
#' ## Import the file using rtracklayer and use the UCSC naming standard
#' library(rtracklayer)
#' surveyed <- import(surveyed_file)
#' seqlevelsStyle(surveyed) <- "UCSC"
#'
#' ## For this example we use the same surveyed file for each sample.
#' surveyed_list <- rep(list(surveyed), 9)
#'
#' ## Calculate the number of observed and expected number of mutations in
#' ## each genomic regions for each sample.
#' distr <- genomic_distribution(vcfs, surveyed_list, regions)
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

genomic_distribution <- function(vcf_list, surveyed_list, region_list) {

  # Check arguments
  if (length(vcf_list) != length(surveyed_list)) {
    stop("vcf_list and surveyed_list must have the same length", call. = FALSE)
  }

  if (is.null(names(region_list))) {
    stop(paste("Please set the names of region_list using:",
      "    names(region_list) <- c(\"regionA\", \"regionB\", ...)",
      sep = "\n"
    ), call. = FALSE)
  }

  # Perform intersect with region over each combi of vcf and region.
  # Map over the region list. Within this loop map over the vcfs.
  df <- purrr::map(as.list(region_list), function(region) {
    purrr::map2(as.list(vcf_list), as.list(surveyed_list), .intersect_with_region, region) %>%
      dplyr::bind_rows(.id = "sample")
  }) %>% dplyr::bind_rows(.id = "region")

  # Region as factor
  # make sure level order is the same as in region_list input (important
  # for plotting later)
  df$region <- factor(df$region, levels = names(region_list))
  return(df)
}
