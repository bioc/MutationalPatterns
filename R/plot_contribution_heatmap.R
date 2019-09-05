#' Plot signature contribution heatmap
#' 
#' Plot relative contribution of signatures in a heatmap
#' 
#' @param contribution Signature contribution matrix
#' @param sig_order Character vector with the desired order of the signature names for plotting. Optional.
#' @param mut_type Character stating which mutation type(s) must be plotted. Values of 'mut_type' must be
#' names of the list 'contribution'
#' @param cluster_samples Hierarchically cluster samples based on eucledian distance. Default = T.
#' @param cluster_mut_type optional signatures subset parameter for clustering
#' @param method The agglomeration method to be used for hierarchical clustering. This should be one of 
#' "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) 
#' or "centroid" (= UPGMC). Default = "complete".
#' @param plot_values Plot relative contribution values in heatmap. Default = F.
#' 
#' @return Heatmap with relative contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom cowplot plot_grid
#'
#' @examples
#' ## See the 'mut_matrix()' example for how we obtained the following
#' ## mutation matrix.
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Extracting signatures can be computationally intensive, so
#' ## we use pre-computed data generated with the following command:
#' 
#' # nmf_res <- extract_signatures(mut_mat, rank = 2)
#'
#' nmf_res <- readRDS(system.file("states/nmf_res_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## Set signature names as row names in the contribution matrix
#' rownames(nmf_res$contribution) = c("Signature A", "Signature B")
#' 
#' ## Define signature order for plotting
#' sig_order = c("Signature B", "Signature A")
#'
#'
#' ## Contribution heatmap with signatures in defined order
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"))
#' 
#' ## Contribution heatmap without sample clustering
#' plot_contribution_heatmap(nmf_res$contribution, 
#'                           sig_order = c("Signature B", "Signature A"), 
#'                           cluster_samples = FALSE, method = "complete")
#'
#' @seealso
#' \code{\link{extract_signatures}},
#' \code{\link{mut_matrix}},
#' \code{\link{plot_contribution}}
#'
#' @export

# plotting function for relative contribution of signatures in heatmap
plot_contribution_heatmap = function(contribution, sig_order, mut_type, cluster_samples = TRUE, cluster_mut_type, method = "complete", plot_values = FALSE)
{
  # check contribution argument
  if(class(contribution) == "list"){combined = F}
  else if (class(contribution) == "matrix")
  {
    warning("Matrix given for 'contribution', treated as combined signatures", call.=T, immediate.=T)
    combined = T
  } else {stop("contribution must be a named list")}
  
  if (!combined)
  {
    # Check mut_type argument
    mut_type = check_mutation_type(mut_type)
    
    if (all(mut_type %in% names(contribution))) {contribution = contribution[mut_type]}
    else {stop(paste("One or more values of 'mut_type' is not found in 'contribution'.", 
                     "Run function without 'mut_type' argument or give values which are in contribution list"))}
    
    if (cluster_samples)
    { 
      if (missing(cluster_mut_type)) {cluster_mut_type = names(contribution)}
      else if (length(cluster_mut_type ) == 1)
      {
        if (cluster_mut_type == "all") {cluster_mut_type = c("snv","dbs","indel")}
      }
      if (!all(cluster_mut_type %in% names(contribution)))
      {stop(paste("One or more values of 'cluster_mut_type' is not found in 'contribution'.", 
                  "Run function without 'cluster_mut_type' argument or give values which are in contribution list"))}
      
      cluster_mutations = c()
      
      for (m in cluster_mut_type)
      {
        if (any(grepl(m, names(contribution)))) { cluster_mutations = c(cluster_mutations, rownames(contribution[[m]])) }
        else { cluster_mutations = cluster_mutations }
      }
      
      if (isEmpty(cluster_mutations))
        warning(paste("Values of 'cluster_mut_type' are not found in names of contribution list or do not match 'mut_type'.",
                      "No clustering on mutation type is done"))
    }  
    
    contribution = do.call(rbind, contribution)
  } 
  
  # check if there are signatures names in the contribution matrix
  if(is.null(row.names(contribution)))
    {stop("contribution must have row.names (signature names)")}
  # if no signature order is provided, use the order as in the input matrix
  if(missing(sig_order))
  {
    sig_order = rownames(contribution)
  }
  # check sig_order argument
  if(class(sig_order) != "character")
    {stop("sig_order must be a character vector")}
  if(length(sig_order) != nrow(contribution))
    {stop("sig_order must have the same length as the number of signatures in the contribution matrix")}
  if(any(is.na(match(sig_order, row.names(contribution)))))
    {stop("sig_order must have the same signature names as in contribution")}
    
  # transpose
  contribution = t(contribution)
  # relative contribution
  contribution_norm = contribution / rowSums(contribution)
  
  # if cluster samples is TRUE, perform clustering
  if (cluster_samples)
  {
    # hiearchically cluster samples based on eucledian distance between relative contribution
    if (combined) {hc.sample = hclust(dist(contribution_norm), method = method)}
    else {hc.sample = hclust(dist(contribution_norm[,match(cluster_mutations, colnames(contribution_norm))]), method = method)}
    # order of samples according to hierarchical clustering
    sample_order = rownames(contribution)[hc.sample$order]
    
    sig_order = c(cluster_mutations, " ", colnames(contribution)[which(!(colnames(contribution) %in% cluster_mutations))])
  } 
  else
  {
    sample_order = rownames(contribution)
  }

  Signature = NULL
  Sample = NULL
  Contribution = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL

  # melt data frame
  contribution_norm.m = melt(contribution_norm)
  # assign variable names
  colnames(contribution_norm.m) = c("Sample", "Signature", "Contribution")
  # change factor levels to the order for plotting
  contribution_norm.m$Sample = factor(contribution_norm.m$Sample, levels = sample_order)
  contribution_norm.m$Signature = factor(contribution_norm.m$Signature, levels = sig_order)
  
  if (cluster_samples)
  {
    if(!isEmpty(cluster_mutations))
      contribution_norm.m = rbind(contribution_norm.m,
                                  data.frame("Sample" = unique(contribution_norm.m$Sample),
                                             "Signature" = " ",
                                             "Contribution" = NA))
  }
    
  # plot heatmap
  heatmap = ggplot(contribution_norm.m, aes(x=Signature, y=Sample, fill=Contribution, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \ncontribution", limits = c(0,1),
                         na.value = "white") +
    scale_x_discrete(breaks=sig_order, labels=sig_order) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x=NULL, y=NULL) 
  # if plot_values is TRUE, add values to heatmap
  if (plot_values)
  {
    heatmap = heatmap + geom_text(aes(label = round(Contribution, 2)), size = 3)
  }
  
  # if cluster_samples is TRUE, make dendrogram
  if (cluster_samples)
  {
    # get dendrogram
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    # plot dendrogram of hierachical clustering
    dendrogram = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      scale_y_reverse(expand = c(0.2, 0)) + 
      theme_dendro()
    # combine plots
    plot_final = cowplot::plot_grid(dendrogram, heatmap, align = 'h', rel_widths = c(0.3, 1))
  } 
  else
  {
    plot_final = heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(contribution_norm.m$Sample))))
  }

  return(plot_final)
}
