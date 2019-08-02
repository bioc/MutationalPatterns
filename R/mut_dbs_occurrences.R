#' Count double base substitutions occurrences
#'  
#' @param type_context result of double base substitutions from type_context function
#' @importFrom S4Vectors isEmpty
#' @noRd
#' @return vector with 78 double base substitutions mutation occurrences

mut_dbs_occurrences = function(type_context)
{
  vector = rep(0,78)
  names(vector) = DBS
  
  # if type_context is empty, return vector with zeroes
  if (isEmpty(type_context))
    return(vector)
  
  # for all mutations in this sample
  for (i in 1:length(type_context))
  {
    # Find mutation type
    type = which(DBS == type_context[i])
    
    # Increase count for mutation type with one
    vector[type] = vector[type] + 1
  }
  
  return(vector)
}
