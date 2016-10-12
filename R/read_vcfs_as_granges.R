#' Read VCF files into a GRangesList
#'
#' This function reads Variant Call Format (VCF) files into a GRanges object
#' and combines them in a GRangesList.  In addition to loading the files, this
#' function applies a seqlevel style to the GRanges objects.  The default
#' seqlevel style is "UCSC".
#'
#' @param vcf_files Character vector of vcf file names
#' @param sample_names Character vector of sample names
#' @param genome A character or Seqinfo object
#' @param style The naming standard to use for the GRanges. (default = "UCSC")
#' @return A GRangesList containing the GRanges obtained from vcf_files
#'
#' @importFrom BiocGenerics match
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb "seqlevelsStyle<-"
#'
#' @examples
#' # The example data set consists of three colon samples, three intestine
#' # samples and three liver samples.  So, to map each file to its appropriate
#' # sample name, we create a vector containing the sample names:
#' sample_names <- c("colon1", "colon2", "colon3",
#'                     "intestine1", "intestine2", "intestine3",
#'                     "liver1", "liver2", "liver3")
#'
#' # We assemble a list of files we want to load.  These files match the
#' # sample names defined above.
#' vcf_files <- list.files(system.file("extdata", 
#'                                     package="MutationalPatterns"),
#'                                     pattern = ".vcf", full.names = TRUE)
#'
#' # This function loads the files as GRanges objects
#' vcfs <- read_vcfs_as_granges(vcf_files, sample_names, genome = "hg19")
#'
#' @export

read_vcfs_as_granges <- function(vcf_files, sample_names, genome="-", 
                                    style="UCSC")
{
    # Check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Provide the same number of sample names as VCF files")

    vcf_list <- GRangesList(lapply (vcf_files, function (file)
    {
        # Use VariantAnnotation's readVcf, but only store the
        # GRanges information in memory.  This speeds up the
        # loading significantly.
        vcf <- rowRanges(readVcf (file, genome))

        # Convert to a single naming standard.
        seqlevelsStyle(vcf) <- style

        # Find and exclude positions with indels or multiple
        # alternative alleles.
        rem <- which(all(!( !is.na(match(vcf$ALT, DNA_BASES)) &
                            !is.na(match(vcf$REF, DNA_BASES)) &
                            (lengths(vcf$ALT) == 1) )))

        if (length(rem) > 0)
        {
            vcf = vcf[-rem]
            warning(length(rem),
                    " position(s) with indels and multiple",
                    " alternative alleles are removed.")
        }

        return(vcf)
    }))

    # Set the provided names for the samples.
    names(vcf_list) <- sample_names

    return(vcf_list)
}

##
## Deprecated variants
##

read_vcf <- function(vcf_files, sample_names, genome="-", style="UCSC")
{
    .Defunct("read_vcfs_as_granges", package="MutationalPatterns",
            msg=paste("This function has been removed.  Use",
                        "'read_vcfs_as_granges' instead.  The new function",
                        "automatically renames the seqlevel style for you,",
                        "so you no longer need to run 'rename_chrom' either."))
}

vcf_to_granges <- function(vcf_files, sample_names, genome="-", style="UCSC")
{
    # Show the same error message as 'read_vcf()'.
    read_vcf()
}

rename_chrom <- function(granges, style = "UCSC")
{
    .Defunct("rename_chrom", package="MutationalPatterns",
            msg = paste("This function has been removed."))
}
