#' Count indel occurrences
#'  
#' @param type_context result of indels from type_context function
#' @param cosmic Boolean stating if indels must be count according
#' to classification in COSMIC database. Default = FALSE
#' 
#' @importFrom S4Vectors isEmpty
#' 
#' @noRd
#' 
#' @return vector with either 30 indel contexts according to own
#' classification or 83 indel contexts according to COSMIC
#' classification

mut_indel_occurrences = function(type_context, cosmic = T)
{
  # if type_context is empty, return vector with zeroes
  if (isEmpty(type_context))
    return(vector)
  
  if (cosmic)
  {
    vector = rep(0,83)
    names(vector) = DBS
    
  } else 
  {
    # for all mutations in this sample
    for (i in 1:length(type_context[[1]]))
    {
      # Find mutation type
      type = which(DBS == type_context[[1]][i])
      
      # Increase count for mutation type with one
      vector[type] = vector[type] + 1
    }
  }
  
 
  
  
  
  return(vector)
}
