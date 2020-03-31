#' Throw an error for when an argument should be a GRanges/GRangesList object.
#'
#' @details 
#' This function is called by other functions, when an argument should be a GRanges/GrangesList object, but isn't.
#' It throws an error showing the actual class of the argument.
#' 
#' @param arg Argument. Any object that should have been a GRanges/GrangesList, but isn't.
#' 
#' @examples
#' \dontrun{
#' a = 1
#' not_gr_or_grl(a)
#' }
#' 

not_gr_or_grl = function(arg){
    arg_name = deparse(substitute(arg))
    arg_class = class(arg)[[1]]
    stop(str_c(arg_name, " should be a CompressedGRangesList or GRanges object, instead it is a ", arg_class, " object.
               Please provide an object of the correct class."))
}