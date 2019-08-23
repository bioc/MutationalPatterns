#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving the nonnegative least-squares
#' constraints problem.
#' 
#' @param mut_matrix Named list of count matrices 
#' @param signatures Named list of signature matrices (number of mutational features
#' for each signature matrix must be the same as in the corresponding count matrix)
#' @param cutoff Numeric value of absolute signature contribution. Signatures 
#' greater than or equal to the given value are returned. Default = 0
#' @param mode Character to select mutation type for which the function has to
#' be executed
#' @param method Character to select the method used to fit the signatures:
#' \itemize{
#'   \item{'least-squares'} {Solve the nonnegative least-squares constraints problem}
#'   \item{'golden-ratio-search'} {Perform gold ratio search algortihm from deconstructSigs package}
#' } 
#' Default is 'least-squares'
#' @param ... Other arguments passed to whichSignatures() when using method = 'golden-ratio-search'
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @examples
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## You can download the signatures from the COSMIC website:
#' # http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#' 
#' ## Match the order to MutationalPatterns standard of mutation matrix
#' order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
#' ## Reorder cancer signatures dataframe
#' cancer_signatures = cancer_signatures[order,]
#' ## Use trinucletiode changes names as row.names
#' ## row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
#' ## Keep only 96 contributions of the signatures in matrix
#' cancer_signatures = as.matrix(cancer_signatures[,4:33])
#' ## Rename signatures to number only
#' colnames(cancer_signatures) = as.character(1:30)
#'
#'
#' ## Perform the fitting
#' fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#' \code{\link{least_square_error_fitting}}
#' \code{\link{golden_ratio_search_fitting}}
#'
#' @export

fit_to_signatures = function(mut_matrix, signatures, mode, cutoff, method = "least-squares", ...)
{
    mode = check_mutation_type(mode)
    
    if (class(signatures) == "matrix")
    { 
      if (class(mut_matrix) == "matrix") { signatures = list("unknown"=signatures) }
      else 
      {
        signatures_list = list()
        for (m in names(mut_matrix))
        {
          if (all(rownames(mut_matrix[[m]]) == rownames(signatures))) 
          { 
            signatures_list[[m]] = signatures
            signatures = signatures_list
            break 
          }
        }
      }
    }
    
    if (class(mut_matrix) == "matrix") { mut_matrix = list("unknown"=mut_matrix) }
    
    mut_matrix = mut_matrix[intersect(mode, names(mut_matrix))]
    
    dots = list(...)
    
    if (method == "least-squares")
    {
      if (missing(cutoff)) { cutoff = 0 }
      res = least_squares_error_fitting(mut_matrix, signatures, cutoff)
    } else if (method == "golden-ratio-search")
    {
      if (!("signature.cutoff" %in% names(dots)) & !missing(cutoff)) 
        { res = golden_ratio_search_fitting(mut_matrix, signatures, signature.cutoff = cutoff, ...) }
      else 
        { res = golden_ratio_search_fitting(mut_matrix, signatures, ...) }
    }
    
    return(res)
}
