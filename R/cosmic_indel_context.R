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
    len = indel$indel_len
    sq = as.character(indel$indel_seq)
    type = indel$indel_type
    n = indel$repeats
    bimh = indel$bimh
    
    if (sq == "A" | sq == "G") 
        sq = as.character(complement(DNAString(sq)))
    
    if (type == "del")
    {
        if (n >= 2)
        {
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
            
            if (n >= 6) { n = "6+" }
            else { n = as.character(n) }
          
            return(sprintf("del.%s.%s", len, n))
        } else {
            if (len == 1)
            {
              return(sprintf("del.1bp.homopol.%s.len.1", sq))
            }
            if (len >= 5) { len = "5+"}
          
            if (bimh >= 1)
            {
                if (bimh >= 5){ bimh = "5+" }
                else { bimh = as.character(bimh) }
                
                return(sprintf("del.mh.len.%s.bimh.%s", len, bimh))
            } else {
                return(sprintf("del.rep.len.%s.rep.1", len))
            }
        }
    } else if (type == "ins") 
    { 
        n = n-1
    
        if (n >= 1)
        {
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
            
            if (n >= 5) { n = "5+" }
            else { n = as.character(n) }
        
            return(sprintf("ins.%s.%s", len, n))
        } else {
            if (len == 1)
            {
                return(sprintf("ins.1bp.homopol.%s.len.0", sq))
            } else if (len >= 5) 
            {
                return("ins.rep.len.5+.rep.0")
            } else {
                return(sprintf("ins.rep.len.%s.rep.0", len))
            }
        }
    }
}
