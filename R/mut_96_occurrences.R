#' Count 96 trinucleotide mutation occurrences
#'  
#' @param type_context result from type_context function
#' @importFrom S4Vectors isEmpty
#' @noRd
#' @return vector with 96 trinucleotide mutation occurrences

mut_96_occurrences = function(type_context)
{
    warning(paste("This function will be deprecated soon,",
                  "use 'mut_occurrences' instead"))

    vector = mut_occurrences(type_context, mode = "snv")
    return(vector)
}
