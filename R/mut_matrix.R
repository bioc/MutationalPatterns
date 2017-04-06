#' Make mutation count matrix of 96 trinucleotides 
#'  
#' @description Make 96 trinucleotide mutation count matrix
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @return 96 mutation count matrix
#' @import GenomicRanges
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#' 
#' ## Exclude mitochondrial and allosomal chromosomes.
#' autosomal <- extractSeqlevelsByGroup(species="Homo_sapiens",
#'                                         style="UCSC",
#'                                         group="auto")
#'
#' vcfs <- lapply(vcfs, function(x) keepSeqlevels(x, autosomal))
#'
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#'
#' ## Et voila.
#' mut_mat
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_matrix = function(vcf_list, ref_genome)
{
    df = data.frame()

    num_cores <- detectCores()
    if (is.na(num_cores)) num_cores = 2

    rows <- mclapply (as.list(vcf_list), function (vcf)
    {
        type_context = type_context(vcf, ref_genome)
        row = mut_96_occurrences(type_context)
        return(row)
    }, mc.cores = (num_cores - 1))

    # Merge the rows into a dataframe.
    for (row in rows)
        df = rbind (df, row)

    names(df) = names(row)
    row.names(df) = names(vcf_list)

    # transpose
    return(t(df))
}
