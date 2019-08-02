#' Compute all pairwise cosine similarities between mutational profiles/signatures
#' 
#' Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices. 
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
#' 
#' @param mut_matrix1 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param mut_matrix2 96 mutation count matrix (dimensions: 96 mutations X m samples)
#' @return Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
#'
#' @examples
#' ## You can download mutational signatures from the COSMIC database:
#' # sp_url = http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                    package="MutationalPatterns"))
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
#' ## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
#' cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}
#' 
#' @export

cos_sim_matrix = function(mut_matrix1, mut_matrix2, mode)
{
  
  if (missing(mode) & class(mut_matrix1) == "matrix" & class(mut_matrix2) == "matrix") { mode = "unknown" }
  else if (missing(mode)) { mode = c("snv","dbs","indel") }
  
  if (class(mut_matrix1) == "matrix" & class(mut_matrix2) == "list")
  { 
    row_list = list()
    for (m in names(mut_matrix2))
    {
      if (all(rownames(mut_matrix2[[m]]) == rownames(mut_matrix1))) 
      { 
        row_list[[m]] = mut_matrix1
        mut_matrix1 = row_list
        break 
      }
    }
  } else if (class(mut_matrix1) == "list" & class(mut_matrix2) == "matrix")
  {
    row_list = list()
    for (m in names(mut_matrix1))
    {
      if (all(rownames(mut_matrix1[[m]]) == rownames(mut_matrix2))) 
      { 
        row_list[[m]] = mut_matrix2
        mut_matrix2 = row_list
        break 
      }
    }
  }
  if (class(mut_matrix1) == "matrix" & class(mut_matrix2) == "matrix") 
  {
    mut_matrix1 = list("unknown"=mut_matrix1)
    mut_matrix2 = list("unknown"=mut_matrix2)
  }
  
  res_matrix_list = list()
  
  names_lists = intersect(names(mut_matrix1), names(mut_matrix2))
  if (all(is.na(names_lists))) { stop("No matching names between 'mut_matrix1' and 'mut_matrix2'")}
  if (any(is.na(names_lists))) { warning("Only cosine similarity matrices for matching names in the mutation matrices")}
  
  for (m in names_lists)
  {
    if (!(m %in% names(mut_matrix2))) { stop("One or more names of 'mut_matrix1' not found in 'mut_matrix2'")}
    
    n_samples1 = ncol(mut_matrix1[[m]])
    n_samples2 = ncol(mut_matrix2[[m]])
    res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
    
    for(s in 1:n_samples1)
    {
      signal1 = mut_matrix1[[m]][,s]
      cos_sim_vector = c()
      for(i in 1:n_samples2)
      {
        signal2 = mut_matrix2[[m]][,i]
        cos_sim_vector[i] = cos_sim(signal1, signal2)
      }
      res_matrix[s,] = cos_sim_vector
    }
    rownames(res_matrix) = colnames(mut_matrix1[[m]])
    colnames(res_matrix) = colnames(mut_matrix2[[m]])
    
    res_matrix_list[[m]] = res_matrix
  }

  if (length(mode) > 1)
  {
    mode = intersect(mode[match(names(mut_matrix1),mode)], mode[match(names(mut_matrix2),mode)])
    mode = mode[which(!is.na(mode))]
    
    if (isEmpty(mode)) { stop("Mode is not found in 'mut_matrix1' or 'mut_matrix2'") }
    res_matrix_list = res_matrix_list[mode]
  }
  else if ( mode == "unknown" ) { res_matrix_list = res_matrix_list[[mode]] }
  else {
    mode = unlist(strsplit(mode,"\\+"))
    mode = intersect(mode[which(mode %in% names(mut_matrix1))], mode[which(mode %in% names(mut_matrix2))])
    mode = mode[which(!is.na(mode))]
    
    if (isEmpty(mode)) { stop("Mode is not found in 'mut_matrix1' or 'mut_matrix2'") }
    res_matrix_list = res_matrix_list[mode]
  }  
  return(res_matrix_list)
}


#' This function has been renamed to 'cos_sim_matrix'.
#'
#' @param mut_matrix 96 mutation count matrix (dimensions: 96 mutations X n samples)
#' @param signatures 96 mutation count matrix (dimensions: 96 mutations X m samples)
#'
#' @return Matrix with pairwise cosine similarities
#'
#' @examples
#' ## You can download mutational signatures from the COSMIC database:
#' # sp_url = http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)
#'
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                    package="MutationalPatterns"))
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
#' ## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
#' cos_sim_matrix(mut_mat, cancer_signatures)
#' 
#' @seealso
#' \code{\link{cos_sim_matrix}}
#' \code{\link{mut_matrix}},
#' \code{\link{fit_to_signatures}},
#' \code{\link{plot_cosine_heatmap}}

explained_by_signatures = function(mut_matrix, signatures)
{
  .Defunct("cos_sim_matrix", package="MutationalPatterns",
           msg=paste("This function has been renamed to",
                     "'cos_sim_matrix'."))
}
