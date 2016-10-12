#' Find transcriptional strand of base substitutions in vcf
#'
#' For the positions that are within gene bodies it is determined whether
#' the "C" or "T" base is on the same strand as the gene definition. (Since
#' by convention we regard base substitutions as C>X or T>X.)
#'
#' Base substitions on the same strand as the gene definitions are considered
#' untranscribed, and on the opposite strand of gene bodies as transcribed,
#' since the gene definitions report the coding or sense strand, which is
#' untranscribed.
#'
#' No strand information "-" is returned for base substitutions outside gene
#' bodies, or base substitutions that overlap with more than one gene body.
#'
#' @param vcf GRanges containing the VCF object
#' @param genes GRanges with gene bodies definitions including strand
#' information
#'
#' @return Character vector with transcriptional strand information with
#' length of vcf: "-" for positions outside gene bodies, "U" for
#' untranscribed/sense/coding strand, "T" for
#' transcribed/anti-sense/non-coding strand.
#'
#' @examples
#' ## For this example we need our variants from the VCF samples, and
#' ## a known genes dataset.  See the 'read_vcfs_as_granges()' example
#' ## for how to load the VCF samples.
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' # Exclude mitochondrial and allosomal chromosomes.
#' autosomal = extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                     style="UCSC",
#'                                     group="auto")
#'
#' vcfs = lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' ## You can obtain the known genes from the UCSC hg19 dataset using
#' ## Bioconductor:
#' # source("https://bioconductor.org/biocLite.R")
#' # biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' # library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#'
#' ## For this example, we preloaded the data for you:
#' genes_hg19 <- readRDS(system.file("states/genes_hg19.rds",
#'                         package="MutationalPatterns"))
#'
#' strand_from_vcf(vcfs[[1]], genes_hg19)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

strand_from_vcf = function(vcf, genes)
{
    # Check consistency of chromosome names.
    if (!(all(seqlevels(vcf) %in% seqlevels(genes))))
        stop(paste( "Chromosome names (seqlevels) of vcf and genes Granges",
                    "object do not match. Use the seqlevelsStyle() function",
                    "to rename chromosome names.") )

    # Determine overlap between vcf positions and genes.
    overlap = findOverlaps(vcf, genes)
    overlap = as.data.frame(as.matrix(overlap))
    colnames(overlap) = c('vcf_id', 'gene_body_id')

    # Remove mutations that overlap with multiple genes and therefore cannot
    # be determined whether they are on transcribed or untranscribed strand
    # duplicated mutations.
    dup_pos = overlap$vcf_id[duplicated(overlap$vcf_id)]

    # Index of duplicated mutations
    dup_idx = which(overlap$vcf_id %in% dup_pos)

    # Remove all duplicated (non-unique mapping) mutations.
    if (length(dup_idx) > 0)
        overlap = overlap[-dup_idx,]

    # Subset of mutations in genes
    vcf_overlap = vcf[overlap$vcf_id]

    # Find reference allele of mutations (and strand of reference genome is
    # reported in vcf file).
    ref = vcf_overlap$REF

    # Find the strand of C or T (since we regard base substitutions as
    # C>X or T>X) which mutations have ref allele C or T.
    i = which(ref == "C" | ref == "T")

    # Store mutation strand info in vector.
    strand_muts = rep(0, nrow(overlap))
    strand_muts[i] = "+"
    strand_muts[-i] = "-"

    # Find strand of gene bodies of overlaps.
    strand_genebodies = as.character(strand(genes)[overlap$gene_body_id])

    # Find if mut and gene_bodies are on the same strand.
    same_strand = (strand_muts  == strand_genebodies)

    # Subset vcf object for both untranscribed and transcribed
    # gene definition represents the untranscribed/sense/coding strand
    # if mutation is on same strand as gene, than its untranscribed.
    U_index = which(same_strand == TRUE)

    # If mutation is on different strand than gene, then its transcribed.
    T_index = which(same_strand == FALSE)
    strand = rep(0, nrow(overlap))
    strand[U_index] = "U"
    strand[T_index] = "T"

    # Make vector with all positions in input vcf for positions that do
    # not overlap with gene bodies, report "-".
    strand2 = rep("-", length(vcf))
    strand2[overlap$vcf_id] = strand
    return(strand2)
}

##
## Deprecated variants
##

get_strand <- function(vcf, genes)
{
    .Defunct("strand_from_vcf", package="MutationalPatterns",
                msg=paste("This function has been removed.  Use",
                            "'strand_from_vcf' instead."))
}
