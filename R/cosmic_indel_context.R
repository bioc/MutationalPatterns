#' Retrieve context of base substitutions
#'
#' A function to extract the bases 3' upstream and 5' downstream of the base
#' substitutions from the reference genome
#' @param vcf A Granges object
#' @param ref_genome Reference genome
#' @param mode A character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' @return Character vector with the context of the base substitutions
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings getSeq
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context <- mut_context(vcfs[[1]], ref_genome, mode)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

cosmic_indel_context = function(indel) 
{
    len = indel[5]
    sq = as.character(indel[7])
    type = indel[6]
    n = indel[8]
    con = indel[9]
    bimh = indel[10]
    
    if (con != "mh")
    {
        if (sq == "A" | sq == "G") 
            sq = as.character(complement(DNAString(sq)))
        
        if (len == 1) 
        { 
          len = paste0("1bp.homopol.",sq,".len") 
        } else if (len >= 5) 
        { 
          len = "rep.len.5+.rep" 
        } else 
        { 
          len = paste0("rep.len.",len,".rep") 
        }
        
        if (type == "del"){
            if (n >= 6) { n = "6+" }
            else { n = as.character(n) }
        } else if (type == "ins")
        {
            if (n >= 5) { n = "5+" }
            else { n = as.character(n) }
        }
        
        final_context = sprintf("%s.%s.%s", type, len, n)
        
    } else {
        if (type == "del")
        {
            if (len >= 5){ len = "5+" }
            else { len = as.character(len) }
          
            if (bimh >= 5){ bimh = "5+" }
            else { bimh = as.character(bimh) }
          
            final_context = sprintf("del.mh.len.%s.bimh.%s", len, bimh)
        } else 
        {
          final_context = NA
        }
    }
    
    return(final_context)
}
