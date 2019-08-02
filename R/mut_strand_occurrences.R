#' Compute strand mutation count vector
#' 
#' Compute strand mutation count vector for all mutation types
#'  
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand, 
#' and "-" for unknown
#' @param mode character stating which type of mutation is to be extracted: 
#' 'snv', 'snv+dbs', 'snv+indel', 'dbs', 'dbs+indel', 'indel' or 'all'
#' 
#' @noRd
#' @return A vector with 192 mutation occurrences and 96 trinucleotides
#' for two strands

mut_strand_occurrences = function(type_context, strand, mode)
{
  # get possible strand values
  values = levels(strand[[mode]])
  
  idx1 = which(strand[[mode]] == values[1])
  idx2 = which(strand[[mode]] == values[2])
  
  # get type context for both vcf subsets
  type_context_1 = lapply(type_context[[mode]], function(x) x[idx1])
  type_context_2 = lapply(type_context[[mode]], function(x) x[idx2])
  
  if (mode == "snv")
  {
    # make 96-trinucleotide count vector per set
    vector1 = mut_96_occurrences(type_context_1)
    vector2 = mut_96_occurrences(type_context_2)
    
    # add names
    names_1 = paste(TRIPLETS_96, values[1], sep = "-")
    names_2 = paste(TRIPLETS_96, values[2], sep = "-")
  } else if (mode == "dbs")
  {
    # make DBS count vector per set
    vector1 = mut_dbs_occurrences(type_context_1[[1]])
    vector2 = mut_dbs_occurrences(type_context_2[[1]])
    
    # add names
    names_1 = paste(DBS, values[1], sep = "-")
    names_2 = paste(DBS, values[2], sep = "-")
  } else { return(NULL) }
  
  # combine vectors in alternating fashion
  vector = c(rbind(vector1, vector2))
  names = c(rbind(names_1, names_2))
  names(vector) = names

  return(vector)
}

##
## Renamed function
##

#'
#' This function has been renamed. Use 'mut_strand_occurrences' instead.
#'
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand, 
#' and "-" for unknown
#'
#' @return A vector with 192 mutation occurrences and 96 trinucleotides
#' for two strands
#' 
#' @seealso 
#' \code{\link{mut_strand_occurrences}}
#' 
#' @export

mut_192_occurrences <- function(type_context, strand)
{
  .Defunct("mut_strand_occurrences", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'mut_strand_occurrences' instead."))
}
