#' Cosine similarity function
#' 
#' Calculate the cosine similarity between two vectors of the same length.
#' The cosine similarity is a value between 0 (distinct) and 1 (identical) 
#' and indicates how much two vectors are alike.
#' 
#' @param x Vector 1 of length n
#' @param y Vector 2 of length n
#' @return Cosine similarity value; a value between 0 and 1
#' 
#' @examples 
#' x = c(1.1, 2.1, 0.2, 0.1, 2.9)
#' y = c(0.9, 1.9, 0.5, 0.4, 3.1)
#' cos_sim(x, y)
#' 
#' @export

cos_sim = function(x, y)
{
  res = x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res = as.numeric(res)
  return(res)
}
