#' Set global context, class and color values for indels
#' 
#' A function to set the global values of the context variable, 
#' the class variable and the color variable of indels. For indels
#' there is no single intuitive and naturally constrained set of 
#' mutation types, so a custom set of mutation types can be used
#'
#' @param indel List of vectors for context, class and color of indels.
#' Color vector must have the same length as the class vector, 
#' because in plotting profiles, each class is represented by one color.\cr\cr
#' It is also possible to give a character for predefined variables:
#' \itemize{
#'   \item{"native"} {Represents a indel context of 3 classes per 
#'   deletion and indel: "repetitive region", "microhomology" and 
#'   "none" of these two. Indels have lengths 1 to 5+}
#'   \item{"cosmic"} {Represents the indel context according to the
#'   COSMIC database}
#' }
#' Default is "cosmic"
#' @return The global variables indel_class, indel_class_header and
#' indel_context will be set with this function. Furthermore a 
#' character is returned, which state the chosen option
#'
#' @examples 
#' ## For custom input of native context
#' indel = list("matrix" = mut_matrix,
#'              "context" =c(
#' paste0('del.rep.len.', 1:5),
#' paste0('ins.rep.len.', 1:5),
#' paste0('del.mh.bimh.', 1:5),
#' paste0('ins.mh.bimh.', 1:5),
#' paste0('del.none.len.', 1:5),
#' paste0('ins.none.len.', 1:5)
#' ),
#'              "class"=c(
#' rep('del.rep', 5), rep('ins.rep', 5),
#' rep('del.mh', 5), rep('ins.mh', 5),
#' rep('del.none', 5), rep('ins.none', 5)
#' ),
#'              "colors"=COLORS_INDEL = c(
#' "#F7BF80", "#ED8212", 
#' "#B5D988", "#31A12C", 
#' "#E44A39", "#B81C20"
#' ))
#' 
#' indel_mutation_type(indel)
#'
#' @export

indel_mutation_type <- function(indel)
{
  if (missing(indel)) { indel = "cosmic" }
  if (class(indel) == "character")
  {
    if (indel == "native"){
      indel_context <<- INDEL
      indel_class <<- INDEL_CLASS
      indel_class_header <<- INDEL_CLASS_HEADER
      indel_colors <<- COLORS_INDEL
    } else if (indel == "cosmic") 
    {
      indel_context <<- INDEL_COSMIC
      indel_class <<- INDEL_COSMIC_CLASS
      indel_class_header <<- INDEL_COSMIC_CLASS_HEADER
      indel_colors <<- COLORS_INDEL_COSMIC
    } else { stop(print("Unknown option for indel mutation types",
                        "Provide either 'native' or 'cosmic' character or",
                        "list with indel classes and contexts")) }
  } else if (class(indel) == "list")
  {
    needed = c("matrix", "class","context","colors")
    if (all(needed %in% names(indel)))
    {
      if ("header" %in% names(indel)) { indel_header <<- indel$header }
      else { indel_header <<- NULL }
      indel_matrix <<- indel$matrix
      indel_context <<- indel$context
      indel_class <<- indel$class
      indel_colors <<- indel$colors
      indel = "custom"
    } else 
    {
      stop(sprintf("Custom indel list misses information about: %s", 
                   needed[which(!(needed %in% names(indel)))] ))
    }
  }
  
  return(indel)
}
