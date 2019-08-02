#' Find optimal nonnegative linear combination of mutation signatures to
#' reconstruct the mutation matrix.
#' 
#' Find the linear combination of mutation signatures that most closely
#' reconstructs the mutation matrix by solving the nonnegative least-squares
#' constraints problem.
#' 
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations
#' X n samples)
#' @param signatures Signature matrix (dimensions: 96 mutations
#' X n signatures)
#' @param cutoff Numeric value of absolute signature contribution. Signatures above
#' value are returned
#'
#' @return Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @importFrom pracma lsqnonneg
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
#'
#' @export

fit_to_signatures = function(mut_matrix, signatures, mode, cutoff = 0)
{
  
    if (missing(mode) & class(mut_matrix) == "matrix") { mode = "unknown" }
    else if (missing(mode)) { mode = c("snv","dbs","indel") }
  
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
    
    contribution = list()
    reconstructed = list()
    
    for (m in names(mut_matrix))
    {
      if (!(m %in% names(signatures))) { stop("One or more names of 'mut_matrix' not found in 'signatures'")}
  
      # make sure dimensions of input matrix are correct
      if (dim(mut_matrix[[m]])[1] != dim(signatures[[m]])[1])
          stop(paste("Mutation matrix and signatures input have",
                     "different number of mutational features"))
  
      n_features = dim(mut_matrix[[m]])[1]
      n_samples = dim(mut_matrix[[m]])[2]
      n_signatures = dim(signatures[[m]])[2]
      lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
      lsq_reconstructed = matrix(NA, nrow=n_features, ncol=n_samples)
  
      # Process each sample
      for (i in 1:ncol(mut_matrix[[m]]))
      {
          y = mut_matrix[[m]][,i]
          lsq = lsqnonneg(signatures[[m]], y)
          lsq_contribution[,i] = lsq$x
          lsq_reconstructed[,i] = signatures[[m]] %*% as.matrix(lsq$x) 
      }
  
      # Add row and col names
      sample_names = colnames(mut_matrix[[m]])
      signature_names = colnames(signatures[[m]])
      mut_type_names = rownames(signatures[[m]])
      
      colnames(lsq_contribution) = sample_names
      rownames(lsq_contribution) = signature_names
      
      lsq_contribution = lsq_contribution[which(rowSums(lsq_contribution)>10),]
  
      colnames(lsq_reconstructed) = sample_names
      rownames(lsq_reconstructed) = mut_type_names
      
      contribution[[m]] = lsq_contribution
      reconstructed[[m]] = lsq_reconstructed
    }
    
    if (mode == "unknown") 
    { 
      contribution = contribution[[mode]]
      reconstructed = reconstructed[[mode]]
    } else {
      mode = unlist(strsplit(mode,"\\+"))
      mode = mode[match(names(contribution),mode)]
      
      if (is.na(mode)) { stop("No modes given are found in 'mut_matrix'") }
      contribution = contribution[mode]
      reconstructed = reconstructed[mode]
    }
    
    
    
    res = list(contribution, reconstructed)
    names(res) = c("contribution", "reconstructed")

    return(res)
}
