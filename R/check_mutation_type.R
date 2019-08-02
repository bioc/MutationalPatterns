#' Check the mutation type argument
#' 
#' A function to check if the mutation type argument is a possible option
#'
#' @param mode
#' @return Mutation types 
#'

check_mutation_type <- function(mode)
{
  if (missing(mode)) {mode = c("snv","dbs","indel")}
  
  # Translate mode to lower case
  mode = tolower(mode)
  
  if (any(mode == "all")) {mode = c("snv","dbs","indel")}
  else if (all(unlist(strsplit(mode, "\\+")) %in% c("snv","dbs","indel"))) {mode = unlist(strsplit(mode, "\\+"))}
  else {stop("One or more mutation types given are unknown")}
  
  return(mode)
}
