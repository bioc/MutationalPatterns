#' Find overlap between mutations and a genomic region
#' 
#' Find the number of mutations that reside in genomic region and take
#' surveyed area of genome into account.
#' 
#' @param vcf CollapsedVCF object with mutations
#' @param surveyed GRanges object with regions of the genome that were surveyed
#' @param region GRanges object with genomic region(s)
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @noRd
#' @return A data.frame containing the overlapping mutations for a
#' genomic region.

intersect_with_region = function(vcf, surveyed, region, mode)
{
    if (mode == "snv") { i = which(nchar(as.character(vcf$REF))==1 & nchar(as.character(unlist(vcf$ALT)))==1) }
    else if (mode == "dbs") { i = which(nchar(as.character(vcf$REF))==2 & nchar(as.character(unlist(vcf$ALT)))==2) }
    else if (mode == "indel") { i = which(nchar(as.character(vcf$REF)) != nchar(as.character(unlist(vcf$ALT)))) }
  
    vcf = vcf[i,]
    
    # Check if chromosome names are the same in the objects
    if (seqlevelsStyle(vcf) != seqlevelsStyle(surveyed))
      stop(paste("The chromosome names (seqlevels) of the VCF and the",
                 "surveyed GRanges object do not match."))
    
    if (seqlevelsStyle(region) != seqlevelsStyle(surveyed))
      stop(paste("The chromosome names (seqlevels) of the surveyed and",
                 "the region GRanges object do not match."))
    
    # Mutations on same chromosomes as in surveyed
    vcf = vcf[which(as.character(seqnames(vcf)) %in% unique(as.character(seqnames(surveyed)))),]
    
    # Number of mutations in vcf file
    n_muts = length(vcf)

    # Number of base pairs that were surveyed
    surveyed_length = sum(as.numeric(width(surveyed)))

    # Intersect genomic region and surveyed region
    surveyed_region = intersect(surveyed, region, ignore.strand = TRUE)
    surveyed_region_length = sum(width(surveyed_region))

    # Find which mutations lie in surveyed genomic region
    overlap = findOverlaps(vcf, surveyed_region)
    muts_in_region = as.data.frame(as.matrix(overlap))$queryHits

    observed = length(muts_in_region)
    prob = n_muts / surveyed_length
    expected = prob * surveyed_region_length

    res = data.frame(n_muts,
                        surveyed_length,
                        prob, surveyed_region_length,
                        expected,
                        observed)
    return(res)
}
