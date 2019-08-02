#' Retrieve mutations from vcf
#' 
#' A function to extract mutations of each position in vcf
#' @param vcf A CollapsedVCF object
#' @param mode A character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' @return List with character vector for each mutation type
#' @import GenomicRanges
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' muts = mutations_from_vcf(vcfs[[1]])
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}}
#'
#' @export

mutations_from_vcf = function(vcf, mode="snv") 
{
    # Translate mode to lower case
    mode = tolower(mode)
    
    if(!(mode %in% c("snv", "snv+dbs","snv+indel","dbs","dbs+indel","indel","all")))
        stop("Mode must be 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'")
    
    ref = as.character(vcf$REF)
    alt = as.character(unlist(vcf$ALT))

    # Allow both uppercase and lowercase column names.
    if (length(ref) == 0)
        ref = as.character(vcf$ref)

    if (length(alt) == 0)
        alt = as.character(vcf$alt)

    # If these columns are still missing, there's nothing we can do.
    if (length(ref) == 0 || length(alt) == 0)
        warning("Some of your data is missing an ALT and/or a REF column.")
    
    muts = list()

    if (mode == "all")
    {
      mode = "snv+dbs+indel"
    }
    
    if (grepl("snv",mode))
    {
      ref_snv = ref[nchar(ref)==1 & nchar(alt)==1]
      alt_snv = alt[nchar(ref)==1 & nchar(alt)==1]
      
      muts = c(muts,list("snv"=paste(ref_snv,alt_snv,sep=">")))
    } 
    if (grepl("dbs",mode))
    {
      ref_dbs = ref[nchar(ref)==2 & nchar(alt)==2]
      alt_dbs = alt[nchar(ref)==2 & nchar(alt)==2]
      
      muts = c(muts,list("dbs"=paste(ref_dbs,alt_dbs,sep=">")))
    }
    if (grepl("indel",mode))
    {
      ref_ind = ref[nchar(ref)!=nchar(alt)]
      alt_ind = alt[nchar(ref)!=nchar(alt)]
      
      muts = c(muts,list("indel"=paste(ref_ind,alt_ind,sep=">")))
    }

    return(muts)
}
