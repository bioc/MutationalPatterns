#' Plot cosine similarity heatmap
#' 
#' Plot pairwise cosine similarities in a heatmap.
#' 
#' 
#' @param cos_sim_matrix Matrix with pairwise cosine similarities.
#'                       Result from \code{\link{cos_sim_matrix}}
#' @param col_order Character vector with the desired order of the columns names for plotting. Optional.
#' @param cluster_rows Hierarchically cluster rows based on eucledian distance. Default = TRUE.
#' @param cluster_cols Hierarchically cluster cols based on eucledian distance. Default = FALSE.
#' @param method The agglomeration method to be used for hierarchical clustering. This should be one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
#' or "centroid" (= UPGMC). Default = "complete".
#' @param plot_values Plot cosine similarity values in heatmap. Default = FALSE.
#'
#' @return Heatmap with cosine similarities
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#' @importFrom magrittr %>% 
#' @examples
#' 
#' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' # You can download the signatures from the COSMIC website:
#' # https://cancer.sanger.ac.uk/cosmic/signature
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/snv_signatures_probabilities.txt",
#'                      package="MutationalPatterns")
#' signatures <- read.table(filename, sep = "\t", header = TRUE)
#' signatures = as.matrix(signatures[,-c(1,2)])
#'
#' ## Calculate the cosine similarity between each signature and each 96 mutational profile
#' cos_matrix = cos_sim_matrix(mut_mat, signatures)
#' 
#' 
#' ## Plot the cosine similarity between each signature and each sample with hierarchical 
#' ## sample clustering and signatures order based on similarity
#' 
#' plot_cosine_heatmap(cos_matrix, cluster_rows = TRUE, cluster_cols = TRUE, method = "complete")
#' 
#' ## You can also plot the similarity of samples with eachother
#' cos_matrix = cos_sim_matrix(mut_mat, mut_mat)
#' plot_cosine_heatmap(cos_matrix, cluster_rows = TRUE, cluster_cols = TRUE, method = "complete")
#' 
#' @seealso
#' \code{\link{mut_matrix}},
#' \code{\link{cos_sim_matrix}},
#' \code{\link{cluster_signatures}} 
#' 
#' @export

plot_cosine_heatmap = function(cos_sim_matrix, col_order = NA, cluster_rows = TRUE, 
                               cluster_cols = FALSE, method = "complete", plot_values = FALSE)
{
  # check explained argument
  if(! inherits(cos_sim_matrix, "matrix"))
  {stop("cos_sim_matrix must be a matrix")}
  # matrix should have row and colnames
  if(length(colnames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing colnames")}
  if(length(rownames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing rownames")}

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Cosine.sim = Signature = Sample = x = y = xend = yend = NULL

    # if cluster samples is TRUE, perform clustering
  if(cluster_rows == TRUE){
    # cluster samples based on eucledian distance between relative contribution
    hc.sample = hclust(dist(cos_sim_matrix), method = method)
    # order samples according to clustering
    sample_order = rownames(cos_sim_matrix)[hc.sample$order]
    
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_rows = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse(expand = c(0.2, 0)) + 
      theme_dendro()
  }
  else{
    sample_order = rownames(cos_sim_matrix)
  }
  
  
  #If cluster_cols is TRUE perform clustering. Else use supplied col_order or
  #the current column order.
  if (!is_na(col_order) & cluster_cols == TRUE){
    stop("col_order can only be provided when cluster_cols is FALSE", call. = F)
  } else if (!is_na(col_order)){
    # check col_order argument
    if(class(col_order) != "character")
    {stop("col_order must be a character vector")}
    if(length(col_order) != ncol(cos_sim_matrix))
    {stop("col_order must have the same length as the number of signatures in the explained matrix")}
  } else if(cluster_cols == TRUE){
    #Cluster cols
    hc.sample2 = cos_sim_matrix %>% 
      t() %>% 
      dist() %>% 
      hclust(method = method)
    col_order = colnames(cos_sim_matrix)[hc.sample2$order]
    
    dhc = as.dendrogram(hc.sample2)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram_cols = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      theme_dendro() +
      scale_y_continuous(expand = c(0.2, 0))
  } else{
    col_order = colnames(cos_sim_matrix)
  }
  

  # melt
  cos_sim_matrix.m = melt(cos_sim_matrix)
  # assign variable names
  colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  
  # change factor levels to the correct order for plotting
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature, levels = col_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
  # plot heatmap
  heatmap = ggplot(cos_sim_matrix.m, aes(x=Signature, y=Sample, fill=Cosine.sim, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x=NULL, y=NULL)
  # if plot_values is TRUE, add values to heatmap
  if (plot_values == TRUE)
  {
    heatmap = heatmap + geom_text(aes(label = round(Cosine.sim, 2)), size = 3)
  }
  
  # Add dendrogram depending on the clustering of the rows and the columns.
  if (cluster_rows == TRUE & cluster_cols == TRUE){
    empty_fig = ggplot() +
      theme_void()
    plot_final = plot_grid(empty_fig, dendrogram_cols, dendrogram_rows, heatmap, 
                           align = "hv", axis = "tblr", rel_widths = c(0.3, 1), rel_heights = c(0.3, 1))
  }
  else if(cluster_rows == TRUE & cluster_cols == FALSE){
    # combine plots
    plot_final = plot_grid(dendrogram_rows, heatmap, align='h', rel_widths=c(0.3,1))
  } else if (cluster_rows == FALSE & cluster_cols == TRUE){
    plot_final = plot_grid(dendrogram_cols, heatmap, align='v', rel_heights =c(0.3,1)) +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  } else{
    plot_final = heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  return(plot_final)
}
