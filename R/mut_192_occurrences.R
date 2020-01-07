#' Compute 192 mutation count vector
#' 
#' Compute 192 mutation count vector, 96 trinucleotide changes X 2 strands
#'  
#' @param type_context result from type_context function
#' @param strand factor with strand information for each
#' position, for example "U" for untranscribed, "T" for transcribed strand, 
#' and "-" for unknown
#' 
#' @noRd
#' @return A vector with 192 mutation occurrences and 96 trinucleotides
#' for two strands

mut_192_occurrences = function(type_context, strand)
{
  
  warning(paste("This function will be deprecated soon,",
                "use 'mut_strand_occurrences' instead"))
  
  vector = mut_strand_occurrences(type_context, strand, mode = "snv")
  return(vector)
}
