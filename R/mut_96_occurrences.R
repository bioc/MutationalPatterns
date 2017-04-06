#' Count 96 trinucleotide mutation occurrences
#'  
#' @param type_context result from type_context function
#' @noRd
#' @return vector with 96 trinucleotide mutation occurrences

mut_96_occurrences = function(type_context)
{
    vector = rep(0,96)
    names(vector) = TRIPLETS_96

    if (is.null(type_context$types) || is.null(type_context$context))
        return(vector)
    
    # for all mutations in this sample
    for (i in 1:length(type_context[[1]]))
    {
        # Find mutation type
        type = which(SUBSTITUTIONS == type_context[[1]][i])

        # Find triplet
        if(type < 4)
            context = which(C_TRIPLETS == type_context[[2]][i])
        else
            context = which(T_TRIPLETS == type_context[[2]][i])

        pos = (type - 1)*16 + context
        vector[pos] = vector[pos] + 1
    }

    return(vector)
}

mut_96_occurences = function (type_context, strand)
{
    .Defunct("mut_96_occurence", package="MutationalPatterns",
            msg=paste("This function has been renamed to",
                        "'mut_96_occurrences'."))
}
