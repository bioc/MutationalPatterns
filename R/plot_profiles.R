#' Plot profiles for mutations found in mutation matrix
#'
#' Plot relative contribution of profiles
#' @param mut_matrix Mutation matrix of single or double substitution and/or indels
#' @param colors Named list with 6 value color vector "snv" for snv, 10 value color vector "dbs" for dbs and
#' ... value color vector "indel" for indels
#' @param ymax Numeric vector for Y axis maximum value, order is snv, dbs, indel. Can be named
#' vector to specify mutation types. When "method == combine", "ymax" is the maximum of the
#' given values. Default = c(0.2, 0.5, 0.5)
#' @param mut_type Character stating which mutation type(s) must be plotted. Values of 'mut_type' must be
#' names of the list 'mut_matrix'
#' @param method A character stating how profiles must be plotted. "split" for seperate plots, 
#' "combine" for all mutations on same row. Default = "split"
#' @param condensed More condensed plotting format. Default = F.
#' @return profile plot of mutations
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## mutation matrix with transcriptional strand information:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'
#' ## Plot the profiles of three samples
#' plot_profiles(mut_mat[,c(1,4,7)])
#'
#' @seealso
#' \code{\link{mut_matrix}}
#'
#' @export

plot_profiles = function(mut_matrix, colors, ymax, mut_type, method = "split", condensed = FALSE)
{
  # Check if mutation matrix is not empty
  if(all(isEmpty(mut_matrix)))
  {
    stop("Provide a named list for 'mut_matrix' with at least one mutation type")
  }
  
  if (class(mut_matrix) == "matrix")
  {
    if (all(rownames(mut_matrix) %in% TRIPLETS_96)){
      mut_matrix = list("snv"=mut_matrix)
    } else if (all(rownames(mut_matrix) %in% DBS))
    {
      mut_matrix = list("dbs"=mut_matrix)
    } else if (all(rownames(mut_matrix) %in% c(TRIPLETS_96,DBS)))
    {
      mut_matrix = list("snv"=mut_matrix[match(TRIPLETS_96,rownames(nmf_res$signatures)),],
                        "dbs"=mut_matrix[match(DBS,rownames(nmf_res$signatures)),])
      method = "combine"
    }
  }
  
  # Check names of list of mutation matrices
  if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
    stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
  }
  
  if (missing(mut_type)) {mut_type = names(mut_matrix)}
  else if (mut_type == "all") {mut_type = c("snv","dbs","indel")}
  else {mut_type = unlist(strsplit(mut_type, "\\+"))}
  if (all(mut_type %in% names(mut_matrix))) {mut_matrix = mut_matrix[mut_type]}
  else {stop("One or more values of 'mut_type' is not found in 'mut_matrix'")}
  
  context = list()
  substitution = list()
  
  # Check color vector length
  # Colors for plotting
  if(missing(colors)){colors=c(list("snv"=COLORS6),list("dbs"=COLORS10))}
  else {colors = list()}
  
  if(length(colors$snv) != 6){stop("Provide snv colors vector with length 6")}
  if(length(colors$dbs) != 10){stop("Provide dbs colors vector with length 10")}
  
  context = c(context, list("snv"=CONTEXTS_96))
  substitution = c(substitution, list("snv"=rep(SUBSTITUTIONS, each=16)))
  
  # Replace mutated base with dot to get context
  substring(context$snv, 2, 2) = "."
  
  context = c(context, list("dbs"=ALT_DBS))
  substitution_dbs <- unlist(lapply(as.list(SUBSTITUTIONS_DBS), function(sub)
  {
    sub <- unlist(strsplit(sub, ">"))[1]
    l <- length(which(startsWith(DBS, sub)))
    return(rep(paste0(sub,">NN"), l))
  }))
  substitution = c(substitution, list("dbs"=substitution_dbs))
  
  df3 = list()
  
  for (m in names(mut_matrix))
  {
    # Relative contribution
    norm_mut_matrix = apply(mut_matrix[[m]], 2, function(x) x / sum(x) )
    
    # Construct dataframe
    df = data.frame(substitution = substitution[[m]], context = context[[m]])
    rownames(norm_mut_matrix) = NULL
    df2 = cbind(df, as.data.frame(norm_mut_matrix))
    df3 = c(df3, list(melt(df2, id.vars = c("substitution", "context"))))
  }
  
  names(df3) <- names(mut_matrix)
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  
  if (missing(ymax)) {ymax = c("snv"=0.2,"dbs"=0.5,"indel"=0.5)}
  if (isEmpty(names(ymax)))
  {
    if (length(ymax) >= 3)
    {
      ymax = ymax[1:3]
      names(ymax) = c("snv","dbs","indel")
    } else {
      ymax = c(ymax, c(0.2,0.5,0.5)[(1+length(ymax)):3])
      names(ymax) = c("snv","dbs","indel")
    }
    warning("No names given for ymax, order used is 'snv', 'dbs', 'indel'", call. = T, immediate.=T)
  } else if (all(!is.na(match(c("snv","dbs","indel"), names(ymax)))))
  {
    ymax = ymax[c("snv","dbs","indel")]
  } else if (any(!(names(ymax) %in% c("","snv","dbs","indel"))))
  {
    stop("Names of ymax are not 'snv', 'dbs' or 'indel'")
  } else if ("" %in% names(ymax))
  {
    stop("No names found for 1 or more values of ymax, give names for all values")
  }
    
  if (method == "split")
  {
    
    plots = list()
    
    for (m in names(mut_matrix))
    {
      if (condensed)
      {
        plot = ggplot(data=df3[[m]], aes(x=context,
                                    y=value,
                                    fill=substitution,
                                    width=1)) +
          geom_bar(stat="identity", colour="black", size=.2) +
          scale_fill_manual(values=colors[[m]]) +
          facet_grid(variable ~ substitution, scales = "free_x") +
          ylab("Relative contribution") +
          coord_cartesian(ylim=c(0,ymax[m])) +
          scale_y_continuous(breaks=seq(0, ymax[m], 0.1)) +
          # no legend
          guides(fill=FALSE) +
          # white background
          theme_bw() +
          # format text
          theme(axis.title.y=element_text(size=12,vjust=1),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=12),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=9),
                strip.text.y=element_text(size=9),
                panel.grid.major.x = element_blank(),
                panel.spacing.x = unit(0, "lines"))
      } else {
        plot = ggplot(data=df3[[m]], aes(x=context,
                                    y=value,
                                    fill=substitution,
                                    width=0.6)) +
          geom_bar(stat="identity", colour="black", size=.2) + 
          scale_fill_manual(values=colors[[m]]) + 
          facet_grid(variable ~ substitution, scales = "free_x") + 
          ylab("Relative contribution") + 
          coord_cartesian(ylim=c(0,ymax[m])) +
          scale_y_continuous(breaks=seq(0, ymax[m], 0.1)) +
          # no legend
          guides(fill=FALSE) +
          # white background
          theme_bw() +
          # format text
          theme(axis.title.y=element_text(size=12,vjust=1),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=12),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=9),
                strip.text.y=element_text(size=9),
                panel.grid.major.x = element_blank())
      }
      
      if (m == "dbs"){
        plot = plot + xlab("alternative")
      }
      
      plots = c(plots, list(plot))
      
    }
    
    plot = plot_grid(plotlist = plots, ncol = 1) 

  } else if (method == "combine")
  {
    df3 = do.call(rbind, df3)
    colors = unname(unlist(colors[mut_type]))
    
    ymax = max(unname(ymax))
    
    if (condensed)
    {
      plot = ggplot(data=df3, aes(x=context,
                                       y=value,
                                       fill=substitution,
                                       width=1)) +
        geom_bar(stat="identity", colour="black", size=.2) +
        scale_fill_manual(values=colors) +
        facet_grid(variable ~ substitution, scales = "free_x") +
        ylab("Relative contribution") +
        coord_cartesian(ylim=c(0,ymax)) +
        scale_y_continuous(breaks=seq(0, ymax, 0.1)) +
        # no legend
        guides(fill=FALSE) +
        # white background
        theme_bw() +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
              axis.text.y=element_text(size=8),
              axis.title.x=element_blank(),
              axis.text.x=element_text(size=5,angle=90,vjust=0.4),
              strip.text.x=element_text(size=9),
              strip.text.y=element_text(size=9),
              panel.grid.major.x = element_blank(),
              panel.spacing.x = unit(0, "lines"))
    } else {
      plot = ggplot(data=df3, aes(x=context,
                                       y=value,
                                       fill=substitution,
                                       width=0.6)) +
        geom_bar(stat="identity", colour="black", size=.2) + 
        scale_fill_manual(values=colors) + 
        facet_grid(variable ~ substitution, scales = "free_x") + 
        ylab("Relative contribution") + 
        coord_cartesian(ylim=c(0,ymax)) +
        scale_y_continuous(breaks=seq(0, ymax, 0.1)) +
        # no legend
        guides(fill=FALSE) +
        # white background
        theme_bw() +
        # format text
        theme(axis.title.y=element_text(size=12,vjust=1),
              axis.text.y=element_text(size=8),
              axis.title.x=element_blank(),
              axis.text.x=element_text(size=5,angle=90,vjust=0.4),
              strip.text.x=element_text(size=9),
              strip.text.y=element_text(size=9),
              panel.grid.major.x = element_blank())
    }
  }
  
  return(plot)
}
