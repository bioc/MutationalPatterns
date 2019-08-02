#' Make mutation count matrix of 96 trinucleotides 
#'  
#' @description Make 96 trinucleotide mutation count matrix
#' @param vcf_list List of collapsed vcf objects
#' @param ref_genome BSGenome reference genome object 
#' @param mode A character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' @param method Character stating how to use the data. method = "split" will give 
#' results for each mutation type seperately, whereas method = "combine" will give 
#' combined signatures. Default is "split"
#' @param num_cores Number of cores used for parallel processing. If no value
#'                  is given, then the number of available cores is autodetected.
#' @return List with 96 mutation count matrix for single base substitutions, 
#' 78 mutation count matrix for double base substitutions and ... mutation count
#' matrix for indels
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
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Construct a mutation matrix from the loaded VCFs in comparison to the
#' ## ref_genome.
#' mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome, mode)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#'
#' @export

mut_matrix = function(vcf_list, ref_genome, mode = "snv", method = "split", num_cores)
{
    # Check value of method
    if (!(method %in% c("split", "combine"))){ stop("Provide the right value of 'method'. Options are 'split' or 'combine'")}
  
    df = list("snv"=data.frame(), "dbs"=data.frame(), "indel"=data.frame())

    if (missing(num_cores))
    {
        num_cores = detectCores()
        if (!(.Platform$OS.type == "windows" || is.na(num_cores)))
            num_cores <- detectCores()
        else
            num_cores = 1
    }
    
    mode = tolower(mode)
    
    rows <- mclapply (as.list(vcf_list), function (vcf)
    {
        row = list()
        
        if (grepl("snv",mode) | mode == "all")
        {
          type_context = type_context(vcf, ref_genome, mode)
          snv = mut_96_occurrences(type_context$snv)
          row = c(row, list("snv"=snv))
        }
        if (grepl("dbs",mode) | mode == "all")
        {
          type_context = mut_type(vcf, mode)
          dbs = mut_dbs_occurrences(type_context$dbs)
          row = c(row, list("dbs"=dbs))
        }
        
        return(row)
    }, mc.cores = num_cores)

    # Merge the rows into a dataframe of each mutation type.
    for (m in names(rows[[1]]))
    {
      for (row in rows)
      {
      
        if (class (row) == "try-error") stop (row)
        df[[m]] = rbind (df[[m]], row[[m]])
        names(df[[m]]) = names(rows[[1]][[m]])
        
      }
      
      names(df[[m]]) = names(rows[[1]][[m]])
      row.names(df[[m]]) = names(vcf_list)
      df[[m]] <- t(df[[m]])
    }
    
    if (method == "split") { return(df[which(!isEmpty(df))]) }
    else if (method == "combine") { return(do.call(rbind, df[which(!isEmpty(df))]))}
}
