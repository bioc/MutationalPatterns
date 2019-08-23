#' Check the mutation type argument
#' 
#' A function to check if the mutation type argument contains real 
#' mutation types which can be processed by this package
#'
#' @param mode Character stating a combination of mutation types. 
#' If no argument is given, then the function will return 'snv'
#' @return Vector of all possible mutation types from 'mode'
#'
#' @examples 
#' ## Check the combination of mutation types
#' mode = "snv+dbs"
#' check_mutation_type(mode)
#'
#' @export

check_mutation_type <- function(mode)
{
  if (missing(mode)) {mode = c("snv")}
  
  # Translate mode to lower case
  mode = tolower(mode)
  
  if (any(mode == "all")) {mode = c("snv","dbs","indel")}
  else if (all(unlist(strsplit(mode, "\\+")) %in% c("snv","dbs","indel"))) {mode = unlist(strsplit(mode, "\\+"))}
  else {stop("One or more mutation types given are unknown")}
  
  return(mode)
}
