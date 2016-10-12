#' Count 96 trinucleotide mutation occurences
#'  
#' @param type_context result from type_context function
#' @return vector with 96 trinucleotide mutation occurences

mut_96_occurences = function(type_context)
{
    vector = rep(0,96)
    names(vector) = TRIPLETS_96

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
