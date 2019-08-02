#' Extract mutational signatures from 96 mutation matrix using NMF
#'
#' Decomposes trinucleotide count matrix into signatures and contribution of
#' those signatures to the spectra of the samples/vcf files.
#'
#' @param mut_matrix Mutation matrix of single or double substitution and/or indels
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200
#' @param method Character stating how to use the data. method = "split" will give 
#' results for each mutation type seperately, whereas method = "combine" will give 
#' combined signatures. Default is "split"
#' @return Named list of mutation matrix, signatures and signature contribution for
#' each mutation type (method = "split") or all mutation types together (method = "combine")
#'
#' @import NMF
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## This function is computational intensive.
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

extract_signatures = function(mut_matrix, rank, nrun = 200, method = "split")
{
    # Check if mutation matrix is not empty
    if(isEmpty(mut_matrix))
    {
      stop("Provide a named list for 'mut_matrix' with at least one mutation type")
    }
    
    if (class(mut_matrix) == "matrix")
    {
      if (all(rownames(mut_matrix) %in% TRIPLETS_96)) { mut_matrix = list("snv"=mut_matrix) }
      else if (all(rownames(mut_matrix) %in% DBS)) { mut_matrix = list("dbs"=mut_matrix) }
      else 
      { 
        warning("Mutation type of 'mut_matrix' is unknown. Treated as combined mutation types")
        method = "combine"
      }
    }
    
    if (method == "split")
    {
      # Check names of list of mutation matrices
      if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
        stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
      }
      
      nmf_res <- lapply(mut_matrix, function(mut)
      {
        # Add a small pseudocount to avoid features with zero counts.
        mut = as.matrix(mut) + 0.0001
        
        # Make sure the rank_range is valid.
        if (!(rank > 0 & rank == round(rank)))
          stop("Rank should be a positive integer")
        
        if (ncol(mut) < max(rank))
          stop(paste( "The rank should be smaller than the number of",
                      "samples in the input matrix.") )
        
        # Calculate NMF
        res = nmf(mut, rank=rank, method="brunet", nrun=nrun, seed=123456)
        
        # Find signatures and contribution of signatures
        signatures = NMF::basis(res)
        contribution = NMF::coef(res)
        
        # Reconstruct mutation matrix
        reconstructed = signatures %*% contribution
        return(list(signatures = signatures,
                    contribution = contribution,
                    reconstructed = reconstructed))
      })
      signatures = list()
      contribution = list()
      reconstructed = list()
      for (m in names(mut_matrix))
      {
        signatures[[m]] = nmf_res[[m]][["signatures"]]
        contribution[[m]] = nmf_res[[m]][["contribution"]]
        reconstructed[[m]] = nmf_res[[m]][["reconstructed"]]
      }
      nmf_res = list(signatures = signatures,
                     contribution = contribution,
                     reconstructed = reconstructed)
    } else if (method == "combine")
    {
      # Check if mut_matrix is a list or already combined to a matrix
      if (class(mut_matrix) == "list") 
      { 
        # Check names of list of mutation matrices
        if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
          stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
        }
        
        mut_matrix = do.call(rbind, mut_matrix) 
      }
      
      # Add a small pseudocount to avoid features with zero counts.
      mut_matrix = mut_matrix + 0.0001
      
      # Make sure the rank_range is valid.
      if (!(rank > 0 & rank == round(rank)))
        stop("Rank should be a positive integer")
      
      if (ncol(mut_matrix) < max(rank))
        stop(paste( "The rank should be smaller than the number of",
                    "samples in the input matrix.") )
      
      # Calculate NMF
      res = nmf(mut_matrix, rank=rank, method="brunet", nrun=nrun, seed=123456)
      
      # Find signatures and contribution of signatures
      signatures = NMF::basis(res)
      contribution = NMF::coef(res)
      
      # Reconstruct mutation matrix
      reconstructed = signatures %*% contribution
      nmf_res = list(signatures = signatures,
                     contribution = contribution,
                     reconstructed = reconstructed)
    }
    return(nmf_res)
}
