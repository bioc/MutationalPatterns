#' Read VCF files into a GRangesList
#'
#' This function reads Variant Call Format (VCF) files into a GRanges object
#' and combines them in a GRangesList.  In addition to loading the files, this
#' function applies the same seqlevel style to the GRanges objects as the
#' reference genome passed in the 'genome' parameter.
#'
#' @param vcf_files Character vector of VCF file names
#' @param sample_names Character vector of sample names
#' @param genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome of your VCFs
#' @param group Selector for a seqlevel group.  All seqlevels outside
#'              of this group will be removed.  Possible values:
#'              * 'all' for all chromosomes;
#'              * 'auto' for autosomal chromosomes;
#'              * 'sex' for sex chromosomes;
#'              * 'auto+sex' for autosomal + sex chromosomes (default);
#'              * 'circular' for circular chromosomes;
#'              * 'none' for no filtering, which results in keeping all
#'                seqlevels from the VCF file.
#' @param check_alleles logical. If TRUE (default) positions with insertions,
#'              deletions and/or multiple alternative alleles are excluded
#'              from the vcf object, since these positions cannot be analysed
#'              with this package.  This setting can be set to FALSE to speed
#'              up processing time only if the input vcf does not contain any
#'              of such positions, as these will cause obscure errors.
#'
#' @return A GRangesList containing the GRanges obtained from 'vcf_files'
#'
#' @importFrom BiocGenerics match
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb "seqlevelsStyle<-"
#' @importFrom GenomeInfoDb "organism"
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom GenomeInfoDb extractSeqlevelsByGroup
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#' @importFrom plyr llply
#'
#' @examples
#' # The example data set consists of three colon samples, three intestine
#' # samples and three liver samples.  So, to map each file to its appropriate
#' # sample name, we create a vector containing the sample names:
#' sample_names <- c ( "colon1", "colon2", "colon3",
#'                     "intestine1", "intestine2", "intestine3",
#'                     "liver1", "liver2", "liver3" )
#'
#' # We assemble a list of files we want to load.  These files match the
#' # sample names defined above.
#' vcf_files <- list.files(system.file("extdata", 
#'                                     package="MutationalPatterns"),
#'                                     pattern = ".vcf", full.names = TRUE)
#'
#' # Get a reference genome BSgenome object.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library("BSgenome")
#' library(ref_genome, character.only = TRUE)
#'
#' # This function loads the files as GRanges objects
#' vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
#'
#' @export

read_vcfs_as_granges <- function(vcf_files, sample_names, genome = "-",
                                 group = "auto+sex", check_alleles = TRUE)
{
    # Check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Please provide the same number of sample names as VCF files")

    # Check whether the user has adapted to the new behavior of the function.
    if (genome == "-")
        stop(paste("Please pass a reference genome string in the 'genome'",
                   "parameter.  This string can be obtained using",
                   "available.genomes() from the BSgenome package."))

    ref_genome <- base::get(genome)
    ref_organism <- GenomeInfoDb::organism(ref_genome)
    ref_style <- seqlevelsStyle(ref_genome)

    # Name the VCF's genome as the name of the genome build instead of
    # the BSgenome package name.
    genome_name <- genome(ref_genome)[[1]]

    # Check the class of the reference genome
    if (!(class(ref_genome) == "BSgenome"))
        stop("Please provide the name of a BSgenome object.")

    # Detect the number of available cores.  Windows does not support forking,
    # only threading, so unfortunately, we have to set it to 1.
    # On confined OS environments, this value can be NA, and in such
    # situations  we need to fallback to 1 core.
    num_cores = detectCores()
    if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
        num_cores <- detectCores()
    else
        num_cores = 1

    # To be able to print warnings from within the mclapply call,
    # we need to explicitly set this option. See:
    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17122
    original_warn_state = getOption("warn")
    options(warn=1)

    # Show the warning once for all VCF files that are loaded with this
    # call to read_vcfs_as_granges.
    if (!check_alleles)
    {
        warning(paste("check_alleles is set to FALSE.  Make sure your input",
                      "VCF does not contain any positions with insertions,",
                      "deletions or multiple alternative alleles, as these",
                      "positions cannot be analysed with MutationalPatterns",
                      "and cause obscure errors."))
    }

    vcf_list <- GRangesList(mclapply (vcf_files, function (file)
    {
        # Use VariantAnnotation's readVcf, but only store the
        # GRanges information in memory.  This speeds up the
        # loading significantly.
        vcf <- rowRanges(readVcf (file, genome_name))

        # Convert to a single naming standard.
        seqlevelsStyle(vcf) <- ref_style[1]

        groups <- c()
        if (group != "none")
        {
            if (group == "auto+sex")
            {
                groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                                    style = ref_style,
                                                    group = "auto"),
                            extractSeqlevelsByGroup(species = ref_organism,
                                                    style = ref_style,
                                                    group = "sex"))

                # In some cases, the seqlevelsStyle returns multiple styles.
                # In this case, we need to do a little more work to extract
                # a vector of seqlevels from it.
                groups_names <- names(groups)
                if (! is.null(groups_names))
                {
                    # The seqlevels in the groups are now duplicated.
                    # The following code deduplicates the list items, so that
                    # creating a data frame will work as expected.
                    unique_names <- unique(groups_names)
                    groups <- llply(unique_names, function(x) groups[groups_names == x])
                    groups <- llply(groups, unlist, recursive = F)

                    # In case there are multiple styles applied, we only use the first.
                    groups <- unique(as.vector(groups[[1]]))
                }
            }
            else
            {
                groups <- extractSeqlevelsByGroup ( species = ref_organism,
                                                   style = ref_style,
                                                   group = group )
                groups <- unique(as.vector(t(groups)))
            }

            # The provided VCF files may not contain all chromosomes that are
            # available in the reference genome.  Therefore, we only take the
            # chromosomes that are actually available in the VCF file,
            # belonging to the filter group.
            groups <- intersect(groups, seqlevels(vcf))

            # We use 'pruning.mode = "tidy"' to minimize the deleterious effect
            # on variants, yet, remove all variants that aren't in the filter
            # group.  By default, keepSeqlevels would produce an error.
            vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
        }

        if (check_alleles)
        {
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
        }

        return(vcf)
    }, mc.cores = num_cores, mc.silent = FALSE))

    # Reset the option.
    options(warn=original_warn_state)
    
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
