% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cos_sim_matrix.R
\name{cos_sim_matrix}
\alias{cos_sim_matrix}
\title{Compute all pairwise cosine similarities between mutational profiles/signatures}
\usage{
cos_sim_matrix(mut_matrix1, mut_matrix2)
}
\arguments{
\item{mut_matrix1}{mutation count matrix (dimensions: a mutation features X n samples)}

\item{mut_matrix2}{96 mutation count matrix (dimensions: a mutation features X m samples)}
}
\value{
Matrix with pairwise cosine similarities (dimensions: n mutational profiles X m mutational profiles)
}
\description{
Computes all pairwise cosine similarities between the mutational profiles provided in the two mutation count matrices.
The cosine similarity is a value between 0 (distinct) and 1 (identical) and indicates how much two vectors are alike.
}
\examples{
## Get signatures
signatures <- get_known_signatures()

## See the 'mut_matrix()' example for how we obtained the mutation matrix:
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))


## Calculate the cosine similarity between each COSMIC signature and each 96 mutational profile
cos_sim_matrix(mut_mat, signatures)
}
\seealso{
\code{\link{mut_matrix}},
\code{\link{fit_to_signatures}},
\code{\link{plot_cosine_heatmap}}
}
