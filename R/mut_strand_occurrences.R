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
  mode = check_mutation_type(mode)
  
  if (all(names(type_context) %in% c("context", "types")))
  {
    context_96 = do.call(rbind, strsplit(TRIPLETS_96,""))[,c(1,3,7)]
    context_96 = paste0(context_96[,1], context_96[,2], context_96[,3])
    if (all(unique(type_context$context) %in% context_96) &
        !is.null(unique(type_context$context)))
    {
      type_context = list("snv"=list("types"=type_context$types,
                                     "context"=type_context$context))
      mode = "snv"
    }
    else if (all(unique(type_context$types) %in% DBS))
    {
      type_context = list("dbs"=list("types"=type_context$types,
                                     "context"=type_context$context))
      mode = "dbs"
    }
    else if (all(unique(type_context$context) %in% indel_context))
    {
      type_context = list("indel"=list("types"=type_context$types,
                                       "context"=type_context$context))
      mode = "indel"
    }
  }
  
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
    vector1 = mut_occurrences(type_context_1, mode = mode)
    vector2 = mut_occurrences(type_context_2, mode = mode)
    
    # add names
    names_1 = paste(TRIPLETS_96, values[1], sep = "-")
    names_2 = paste(TRIPLETS_96, values[2], sep = "-")
  } else if (mode == "dbs")
  {
    # make DBS count vector per set
    vector1 = mut_occurrences(type_context_1, mode = mode)
    vector2 = mut_occurrences(type_context_2, mode = mode)
    
    # add names
    names_1 = paste(DBS, values[1], sep = "-")
    names_2 = paste(DBS, values[2], sep = "-")
  } else if (mode == "indel")
  { 
    # make indel count vector per set
    vector1 = mut_occurrences(type_context_1, mode = mode, indel = indel_name)
    vector2 = mut_occurrences(type_context_2, mode = mode, indel = indel_name)
    
    # add names
    names_1 = paste(indel_context, values[1], sep = "-")
    names_2 = paste(indel_context, values[2], sep = "-")
  }
  
  # combine vectors in alternating fashion
  if (length(vector1) == 0)
  {
    vector1 = rep(0, length(names_1))
    names(vector1) = names_1
  } 
  if (length(vector2) == 0)
  {
    vector2 = rep(0, length(names_2))
    names(vector2) = names_2
  }
  
  vector = c(rbind(vector1, vector2))
  names = c(rbind(names_1, names_2))
  names(vector) = names

  return(vector)
}