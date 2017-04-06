#' Count 96 trinucleotide mutation occurrences with transcriptional strand
#' information
#'  
#' @param type_context result from type_context function
#' @param strand character vector with strand information for each
#' position, "U" for untranscribed, "T" for transcribed, "-" for
#' unknown or positions outside gene bodies
#' @noRd
#' @return A vector with 192 mutation occurrences and 96 trinucleotides
#' for both transcribed and untranscribed strand

mut_192_occurrences = function(type_context, strand)
{
    U_idx = which(strand == "U")
    T_idx = which(strand == "T")

    # get type context for both vcf subsets
    type_context_U = lapply(type_context, function(x) x[U_idx])
    type_context_T = lapply(type_context, function(x) x[T_idx])

    # make 96-trinucleotide count vector per set
    U_vector = mut_96_occurrences(type_context_U)
    T_vector = mut_96_occurrences(type_context_T)

    # add names
    names_U = paste(TRIPLETS_96, "-u", sep = "")
    names_T = paste(TRIPLETS_96, "-t", sep = "")

    # combine vectors in alternating fashion
    vector = c(rbind(U_vector, T_vector))
    names = c(rbind(names_U, names_T))
    names(vector) = names

    return(vector)
}

mut_192_occurences = function (type_context, strand)
{
    .Defunct("mut_192_occurence", package="MutationalPatterns",
            msg=paste("This function has been renamed to",
                        "'mut_192_occurrences'."))
}
