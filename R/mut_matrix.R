#' Make mutation count matrix of 96 trinucleotides 
#'  
#' @description Make 96 trinucleotide mutation count matrix
#' @param grl GRangesList or GRanges object.
#' @param ref_genome BSGenome reference genome object
#' @param vcf_list Deprecated argument. Replaced with grl
#' @return 96 mutation count matrix
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(grl = grl, ref_genome = ref_genome)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export
mut_matrix = function (grl, ref_genome, vcf_list = NA) {
    
    if (!missing("vcf_list")){
        warning("vcf_list is deprecated, use grl instead. 
              The parameter grl is set equal to the parameter vcf_list.")
        grl <- vcf_list
    }
    
    #Convert list to grl if necessary
    if (inherits(grl, "list")){
        grl = GenomicRanges::GRangesList(grl)
    }
    if (inherits(grl, "CompressedGRangesList")){
        gr_sizes = S4Vectors::elementNROWS(grl)
        gr = unlist(grl)
    } else if (inherits(grl, "GRanges")){
        gr = grl
        gr_sizes = length(gr)
        names(gr_sizes) = "My_sample"
    } else{
        not_gr_or_grl(grl)
    }
    type_context = type_context(gr, ref_genome)
    mut_mat = mut_96_occurrences(type_context, gr_sizes)
    return(mut_mat)
}
