#' Plot 192 trinucleotide profile
#'
#' Plot relative contribution of 192 trinucleotides      
#' @param mut_matrix 192 trinucleotide profile matrix
#' @param ymax Y axis maximum value, default = 0.2
#' @param colors 6 value color vector
#' @param condensed More condensed plotting format. Default = F.
#' @return 192 trinucleotide profile plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## mutation matrix with transcriptional strand information:
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Extract the signatures.
#' ## This is a computationally intensive task, so we load a precomputed
#' ## version instead.
#' # nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
#' nmf_res_strand <- readRDS(system.file("states/nmf_res_strand_data.rds",
#'                                     package="MutationalPatterns"))
#'
#' ## Optionally, provide signature names
#' colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")
#'
#' ## Generate the plot
#' plot_192_profile(nmf_res_strand$signatures)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{extract_signatures}}
#'
#' @export

plot_strand_profiles = function(mut_matrix, colors, ymax, mut_type, mode, method="split", condensed = FALSE)
{
    # Check if mutation matrix is not empty
    if(all(isEmpty(mut_matrix)))
    {
      stop("Provide a named list for 'mut_matrix' with at least one mutation type")
    }
  
    if (!(mode %in% c("transcription", "replication"))) { stop("No option for 'mode' is given. Specify strand bias: 'transcription' or 'replication'")  }
    
    if (class(mut_matrix) == "matrix")
    {
      if (mode == "transcription")
      {
        if (all(rownames(mut_matrix) %in% TRIPLETS_192_trans)){ mut_matrix = list("snv"=mut_matrix) }
        else if (all(rownames(mut_matrix) %in% DBS_trans)) { mut_matrix = list("dbs"=mut_matrix) } 
        else if (all(rownames(mut_matrix) %in% c(TRIPLETS_192_trans,DBS_trans)))
        {
          mut_matrix = list("snv"=mut_matrix[match(TRIPLETS_192_trans,rownames(nmf_res$signatures)),],
                            "dbs"=mut_matrix[match(DBS_trans,rownames(nmf_res$signatures)),])
          method = "combine"
          mut_type = "all"
        }
      } else {
      
        if (all(rownames(mut_matrix) %in% TRIPLETS_192_rep)) { mut_matrix = list("snv"=mut_matrix) }
        else if (all(rownames(mut_matrix) %in% DBS_rep)) { mut_matrix = list("dbs"=mut_matrix) }
        else if (all(rownames(mut_matrix) %in% c(TRIPLETS_192_rep,DBS_rep)))
        {
          mut_matrix = list("snv"=mut_matrix[match(TRIPLETS_192_rep,rownames(nmf_res$signatures)),],
                            "dbs"=mut_matrix[match(DBS_rep,rownames(nmf_res$signatures)),])
          method = "combine"
          mut_type = "all"
        }
      }
    }
    
    # Check names of list of mutation matrices
    if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
      stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
    }
  
    mut_type = check_mutation_type(mut_type)
    
    mut_matrix = mut_matrix[intersect(names(mut_matrix), mut_type)]

    # Check color vector length
    # Colors for plotting
    if(missing(colors)){colors=list("snv"=COLORS6, "dbs"=COLORS10)}
    
    if(length(colors$snv) != 6){stop("Provide colors vector for single base substitutions with length 6")}
    if(length(colors$dbs) != 10){stop("Provide colors vector for double base substitutions with length 10")}
    
    context = list()
    substitution = list()
    strand = list()
    
    context = c(context, list("snv"=rep(CONTEXTS_96, each=2)))
    substitution = c(substitution, list("snv"=rep(SUBSTITUTIONS, each=32)))
    
    # Replace mutated base with dot to get context
    substring(context$snv, 2, 2) = "."
    
    context = c(context, list("dbs"=rep(ALT_DBS, each=2)))
    substitution_dbs <- unlist(lapply(as.list(SUBSTITUTIONS_DBS), function(sub)
    {
      sub <- unlist(strsplit(sub, ">"))[1]
      l <- length(which(startsWith(DBS, sub)))
      return(rep(paste0(sub,">NN"), l))
    }))
    substitution = c(substitution, list("dbs"=rep(substitution_dbs, each=2)))
    
    # These variables will be available at run-time, but not at compile-time.
    # To avoid compiling trouble, we initialize them to NULL.
    value = NULL
    
    if (missing(ymax)) {ymax = c("snv"=NA,"dbs"=NA,"indel"=NA)}
    else if (isEmpty(names(ymax)))
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
      df3 = list()
      
      # get strand from rownames of mut_matrix
      for (m in names(mut_matrix))
      {
        strand[[m]] = sapply(rownames(mut_matrix[[m]]), function(x) strsplit(x, "-")[[1]][2])
        
        norm_mut_matrix = apply(mut_matrix[[m]], 2, function(x) x / sum(x))
        norm_mut_matrix_sum = matrix(nrow=nrow(norm_mut_matrix)/2, ncol=ncol(norm_mut_matrix))
        for (i in 1:nrow(norm_mut_matrix_sum)) 
          norm_mut_matrix_sum[i,] = colSums(norm_mut_matrix[(2*i-1):(2*i),])
        
        if (is.na(ymax[m])) { ymax[m] = max(norm_mut_matrix_sum) }
        
        # Construct dataframe
        df = data.frame(substitution = substitution[[m]],
                        context = context[[m]],
                        strand = strand[[m]])
        
        rownames(norm_mut_matrix) = NULL
        
        df2 = cbind(df, as.data.frame(norm_mut_matrix))
        df3[[m]] = melt(df2, id.vars = c("substitution", "context", "strand"))
      }
      
      plots = list()
      for (m in names(mut_matrix))
      {
        if (ymax[m] <= 0.1) { scale_y = 0.01 }
        else { scale_y = 0.1 }
        
        if (condensed)
        {
          plot = ggplot(data=df3[[m]], aes(x=context,
                                      y=value,
                                      fill=substitution,
                                      width=1,
                                      alpha=strand)) +
            geom_bar(stat="identity", colour="black", size=.2) +
            scale_fill_manual(values=colors[[m]]) +
            facet_grid(variable ~ substitution, scales = "free_x") +
            ylab("Relative contribution") +
            coord_cartesian(ylim=c(0,ymax[m])) +
            scale_y_continuous(breaks=seq(0, ymax[m], scale_y)) +
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
                                    width=0.6,
                                    alpha=strand)) +
            geom_bar(stat="identity", colour="black", size=.2) + 
            scale_fill_manual(values=colors[[m]]) + 
            facet_grid(variable ~ substitution, scales = "free_x") + 
            ylab("Relative contribution") + 
            coord_cartesian(ylim=c(0,ymax[m])) +
            scale_y_continuous(breaks=seq(0, ymax[m], scale_y)) +
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
        
        if (m == "dbs") { plot = plot + xlab('alternative') }
        plots[[m]] = plot
      }
      
      plot = plot_grid(plotlist = plots, ncol = 1)
      return(plot)
    } else 
    {
      colors = unname(unlist(colors[names(mut_matrix)]))
      mut_matrix = do.call(rbind, mut_matrix)
      strand = sapply(rownames(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
      
      norm_mut_matrix = apply(mut_matrix, 2, function(x) x / sum(x))
      norm_mut_matrix_sum = matrix(nrow = nrow(norm_mut_matrix)/2, ncol = ncol(norm_mut_matrix))
      for (i in 1:nrow(norm_mut_matrix_sum))
        norm_mut_matrix_sum[i,] = colSums(norm_mut_matrix[(2*i-1):(2*i),])
      
      if (all(is.na(ymax))) { ymax = max(norm_mut_matrix_sum) }
      else { ymax = max(ymax) }
      
      if (ymax <= 0.1) { scale_y = 0.01 }
      else { scale_y = 0.1 }
      
      # Construct dataframe
      df = data.frame(substitution = unname(unlist(substitution)),
                      context = unname(unlist(context)),
                      strand = strand)
      
      rownames(norm_mut_matrix) = NULL
      
      df2 = cbind(df, as.data.frame(norm_mut_matrix))
      df3 = melt(df2, id.vars = c("substitution", "context", "strand"))
      
      df3$substitution = factor(df3$substitution, levels = unique(unname(unlist(substitution))))
      df3$context = as.character(df3$context)

      if (condensed)
      {
        plot = ggplot(data=df3, aes(x=context,
                                         y=value,
                                         fill=substitution,
                                         width=1,
                                         alpha=strand)) +
          geom_bar(stat="identity", colour="black", size=.2) +
          scale_fill_manual(values=colors) +
          facet_grid(variable ~ substitution, scales = "free_x") +
          ylab("Relative contribution") +
          coord_cartesian(ylim=c(0,ymax)) +
          scale_y_continuous(breaks=seq(0, ymax, scale_y)) +
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
        plot = ggplot(data=df3, aes(x=context,
                                         y=value,
                                         fill=substitution,
                                         width=0.6,
                                         alpha=strand)) +
          geom_bar(stat="identity", colour="black", size=.2) + 
          scale_fill_manual(values=colors) + 
          facet_grid(variable ~ substitution, scales = "free_x") + 
          ylab("Relative contribution") + 
          coord_cartesian(ylim=c(0,ymax)) +
          scale_y_continuous(breaks=seq(0, ymax, scale_y)) +
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
      return(plot)
    }
}
