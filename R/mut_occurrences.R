#' Count occurrences of mutations 
#' 
#' A function to count the occurrences of mutations for all mutation types
#'  
#' @param type_context result of mutations from type_context function
#' @param type (Optional) A character vector stating which type of mutation is to be extracted: 
#' 'snv', 'dbs' and/or 'indel'. All mutation types can also be chosen by 'type = all'.\cr
#' Default is 'snv'
#' @param indel (Optional) A character stating which indel context database to choose:
#' 'predefined' or 'cosmic'. Is used as argument for extract_indels()
#' @importFrom S4Vectors isEmpty
#' @noRd
#' @return List of vector with mutation occurrences for all mutation types asked for

mut_occurrences = function(type_context, type, indel)
{
  # Check mutation type
  type = check_mutation_type(type)
  
  # If type_context only exists of one mutation type, try to find
  # which one
  if (all(names(type_context) %in% c("context", "types")))
  {
    context_96 = do.call(rbind, strsplit(TRIPLETS_96,""))[,c(1,3,7)]
    context_96 = paste0(context_96[,1], context_96[,2], context_96[,3])
    if (all(unique(type_context$context) %in% context_96) &
        !is.null(unique(type_context$context)))
    {
      type_context = list("snv"=list("types"=type_context$types,
                                     "context"=type_context$context))
      type = "snv"
    }
    else if (all(unique(type_context$types) %in% DBS))
    {
      type_context = list("dbs"=list("types"=type_context$types,
                                     "context"=type_context$context))
      type = "dbs"
    }
    else if (all(unique(type_context$context) %in% indel_context))
    {
      type_context = list("indel"=list("types"=type_context$types,
                                     "context"=type_context$context))
      type = "indel"
    }
  }
  count = list()
 
  for (m in type)
  {
    # if type_context is empty, return vector with zeroes for mutation type
    if (all(isEmpty(type_context[[m]])))
    {
      count[[m]] = NULL
      next
    }
    
    # For each mutation type, count the context occurrences
    if (m == "snv")
    {
      count[[m]] = rep(0,96)
      names(count[[m]]) = TRIPLETS_96
      context = sprintf("%s[%s]%s", 
                        substr(type_context[[m]][["context"]], 1, 1),
                        type_context[[m]][["types"]],
                        substr(type_context[[m]][["context"]], 3, 3))
      context = table(context)
    } else if (m == "dbs") 
    {
      count[[m]] = rep(0,78)
      names(count[[m]]) = DBS
      context = table(type_context[[m]][["types"]])
    } else if (m == "indel")
    {
      if (missing(indel)) { indel = "cosmic" }
      if (indel == "cosmic")
      { 
        count[[m]] = rep(0,83) 
        names(count[[m]]) = indel_context
      } else if (indel == "native")
      { 
        count[[m]]=rep(0,30)
        names(count[[m]]) = indel_context
      }
      
      context = table(type_context[[m]][["context"]])
    }
    
    # Give counts to all mutations found. 
    # Some contexts can have a count of 0
    count[[m]][names(context)] = context
  }
  
  if (length(count) == 1){
    count = unlist(unname(count))
  }
  
  return(count)
}
