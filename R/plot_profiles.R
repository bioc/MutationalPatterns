#' Plot profiles for mutations found in mutation matrix
#'
#' Plot relative contribution of mutational profiles
#' @param mut_matrix Names list of mutation count matrices for mutation types. 
#' @param colors Named list with 6 value color vector "snv" for snv, 10 value color vector "dbs" for dbs.
#' For indels give same number of colors as there are classes. \cr
#' 'native' contains 6 classes and 'cosmic' contains 16 classes
#' @param ymax Numeric vector for Y axis maximum value, order is snv, dbs, indel. Can be named
#' vector to specify mutation types. When "method = combine", "ymax" is the maximum of the
#' given values.\cr
#' Set "ymax = maximum" to use the maximal signature contribution for each mutation type as maximum
#' for the Y axis
#' As default, "ymax" is c("snv"=0.2, "dbs"=0.5, "indel"=0.5)
#' @param mut_type Character stating which mutation type(s) must be plotted. Values of "mut_type" must be
#' names of the list "mut_matrix"
#' @param method A character stating how profiles must be plotted. "split" for seperate plots, 
#' "combine" for all mutations on same row. \cr 
#' Default = "split"
#' @param condensed More condensed plotting format. \cr 
#' Default = F.
#' @return profile plot of mutations
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics cbind
#' @importFrom grid unit.c
#' @importFrom gtable gtable_add_rows
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
  
  mut_type = check_mutation_type(mut_type)
  if ("indel" %in% mut_type)
  {
    if (!exists("indel_context")) { stop("Run 'indel_mutation_type()' to set global variables for indels")}
    else if (indel_context[1] == "del.1bp.homopol.C.len.1") { indel = "cosmic" }
    else if (indel_context[1] == "del.rep.len.1") { indel = "native" }
  }
  
  # Check mutation matrix when type is not given
  
  if (class(mut_matrix) == "matrix")
  {
    if (all(rownames(mut_matrix) %in% TRIPLETS_96)){
      mut_matrix = list("snv"=mut_matrix)
    } else if (all(rownames(mut_matrix) %in% DBS))
    {
      mut_matrix = list("dbs"=mut_matrix)
    } else if (all(rownames(mut_matrix) %in% indel_context))
    {
      mut_matrix = list("indel"=mut_matrix)
    } else if (all(rownames(mut_matrix) %in% c(TRIPLETS_96,DBS)))
    {
      mut_matrix = list("snv"=mut_matrix[rownames(mut_matrix) %in% TRIPLETS_96,],
                        "dbs"=mut_matrix[rownames(mut_matrix) %in% DBS,])
      method = "combine"
    } else if (all(rownames(mut_matrix) %in% c(TRIPLETS_96,indel_context)))
    {
      mut_matrix = list("snv"=mut_matrix[rownames(mut_matrix) %in% TRIPLETS_96,],
                        "indel"=mut_matrix[rownames(mut_matrix) %in% indel_context,])
      method = "combine"
    } else if (all(rownames(mut_matrix) %in% c(DBS, indel_context)))
    {
      mut_matrix = list("dbs"=mut_matrix[rownames(mut_matrix) %in% DBS,],
                        "indel"=mut_matrix[rownames(mut_matrix) %in% indel_context,])
      method = "combine"
    } else if (all(rownames(mut_matrix) %in% c(TRIPLETS_96, DBS, indel_context)))
    {
      mut_matrix = list("snv"=mut_matrix[rownames(mut_matrix) %in% TRIPLETS_96,],
                        "dbs"=mut_matrix[rownames(mut_matrix) %in% DBS,],
                        "indel"=mut_matrix[rownames(mut_matrix) %in% indel_context,])
      method = "combine"
    } else {
      stop("Mutation matrix is not a list and mutation types could not be derived")
    }
  }
  
  # Check names of list of mutation matrices
  if (isEmpty(names(mut_matrix)) | any(!(names(mut_matrix) %in% c("snv","dbs","indel")))){
    stop("Provide the right names to the list of mutation matrices. Options are 'snv', 'dbs' and 'indel'")
  }
  
  # Select mutation types which are both in mut_matrix and in mut_type
  mut_type = intersect(names(mut_matrix), mut_type)
  mut_matrix = mut_matrix[intersect(names(mut_matrix), mut_type)]
  
  # Check color vector length
  # Colors for plotting
  if(missing(colors)) { colors=c(list("snv"=COLORS6),list("dbs"=COLORS10),list("indel"=indel_colors)) }
  else {colors = list()}
  
  indel_color_number = 1
  for (i in 2:length(indel_class))
  {
    if (indel_class[i-1] != indel_class[i]) { indel_color_number = indel_color_number + 1 }
  }
  
  if(length(colors$snv) != 6){stop("Provide snv colors vector with length 6")}
  if(length(colors$dbs) != 10){stop("Provide dbs colors vector with length 10")}
  if(length(unique(colors$indel)) != indel_color_number){stop("Provide indel colors vector with length same number of classes")}
  
  # Create context and classes lists
  
  context = list()
  substitution = list()
  
  # Replace mutated base with dot to get context
  substr(CONTEXTS_96, 2, 2) = "."
  context[["snv"]] = CONTEXTS_96
  substitution[["snv"]] = rep(SUBSTITUTIONS, each=16)
  
  # Replace alt with NN in dbs classes
  context[["dbs"]] = ALT_DBS 
  substitution[["dbs"]] <- unlist(lapply(as.list(SUBSTITUTIONS_DBS), function(sub)
  {
    sub <- unlist(strsplit(sub, ">"))[1]
    l <- length(which(startsWith(DBS, sub)))
    return(rep(paste0(sub,">NN"), l))
  }))
  
  if (indel == "cosmic")
  {
    context[["indel"]] = do.call(rbind, strsplit(indel_context, "\\."))[,lengths(strsplit(indel_context, "\\."))[1]]
    df = do.call(rbind, strsplit(indel_context, "\\."))[,c(1:4)]
    indel_class_value = paste(df[,1],df[,2],df[,3],df[,4], sep=".")
    substitution[["indel"]] = indel_class_value
    labels = c("C","T","C","T","2","3","4","5+","2","3","4","5+","2","3","4","5+")
    names(labels) = unique(substitution[["indel"]])
  } else if (indel == "native") {
    context[["indel"]] = do.call(rbind, strsplit(indel_context, "\\."))[,lengths(strsplit(indel_context, "\\."))[1]]
    substitution[["indel"]] = indel_class
  }
  
  context = context[mut_type]
  substitution = substitution[mut_type]
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  
  if (missing(ymax)) {ymax = c("snv"=0.2,"dbs"=0.5,"indel"=0.5)}
  else if (ymax == "maximum") { ymax = c("snv"=NA, "dbs"=NA, "indel"=NA)}
  if (isEmpty(names(ymax)))
  {
    if (length(ymax) >= 3)
    {
      ymax = ymax[1:3]
      names(ymax) = c("snv","dbs","indel")
    } else {
      ymax = c(ymax, c(NA,NA,NA)[(1+length(ymax)):3])
      names(ymax) = c("snv","dbs","indel")
    }
    warning(paste("No names given for ymax, order used is 'snv', 'dbs', 'indel'",
                  "and default ymax is calculated for missing names.",
                  "When 'method'='combine', first value of ymax is used"), call. = T, immediate.=T)
    
  } 
  else if (all(names(mut_matrix) %in% names(ymax))) { ymax = ymax[names(mut_matrix)] } 
  else { stop("None of unknown names found for 1 or more values of ymax, give names of mutation types in 'mut_matrix'") }
    
  if (method == "split")
  {
    
    df3 = list()
    scale_y = c()
    
    for (m in names(mut_matrix))
    {
      # Relative contribution
      norm_mut_matrix = apply(mut_matrix[[m]], 2, function(x) x / sum(x) )
      
      if (is.na(ymax[m])) { ymax[m] = max(norm_mut_matrix) }
      if (ymax[m] <= 0.1) { scale_y[m] = 0.02 }
      else { scale_y[m] = 0.1 }
      
      # Construct dataframe
      df = data.frame(substitution = substitution[[m]], context = context[[m]])
      rownames(norm_mut_matrix) = NULL
      if (m == "indel" & !isEmpty(indel_class_header))
      {
        df2 = cbind(df, "header"=indel_class_header, as.data.frame(norm_mut_matrix))
        df3[[m]] = melt(df2, id.vars = c("header","substitution", "context"))
        df3[[m]]$header = factor(df3[[m]]$header, levels = unique(indel_class_header))
      }
      else 
      {
        df2 = cbind(df, as.data.frame(norm_mut_matrix))
        df3[[m]] = melt(df2, id.vars = c("substitution", "context"))
      }
      df3[[m]]$substitution = factor(df3[[m]]$substitution, levels = unique(df3[[m]]$substitution))
    }
    
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
          ylab("Relative contribution") +
          coord_cartesian(ylim=c(0,ymax[m])) +
          scale_y_continuous(breaks=seq(0, ymax[m], scale_y[m])) +
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
          scale_y_continuous(breaks=seq(0, ymax[m], scale_y[m])) +
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
      
      if (m == "snv"){
        plot = plot + facet_grid(variable ~ substitution, scales = "free_x")
      } else if (m == "dbs"){
        plot = plot + facet_grid(variable ~ substitution, scales = "free_x") + xlab("alternative")
      } else if (m == "indel" & !isEmpty(indel_class_header)){
        plot = plot + 
          facet_nested(variable ~ header + substitution, 
                     scales = "free_x", 
                     labeller = labeller(substitution=labels)) + 
          xlab("length")
      } else {
        plot = plot + facet_grid(variable ~ substitution, scales = "free_x") + xlab("length")
      }
      
      plots = c(plots, list(plot))
      
    }
    
    plot = plot_grid(plotlist = plots, ncol = 1, align = "v", axis = "lr") 

  } else if (method == "combine")
  {
    mut_matrix = do.call(rbind, mut_matrix)
    
    # Relative contribution
    norm_mut_matrix = apply(mut_matrix, 2, function(x) x / sum(x) )
    
    if (is.na(ymax[1])) { ymax = max(norm_mut_matrix) }
    else { ymax = ymax[1] }
    if (ymax <= 0.1) { scale_y = 0.02 }
    else { scale_y = 0.1 }
    
    # Construct dataframe
    df = data.frame(substitution = unname(unlist(substitution)), 
                    context = unname(unlist(context)))
    rownames(norm_mut_matrix) = NULL
    df2 = cbind(df, as.data.frame(norm_mut_matrix))
    df3 = melt(df2, id.vars = c("substitution", "context"))
    
    df3$substitution = factor(df3$substitution, levels = unique(unname(unlist(substitution))))
    
    colors = unname(unlist(colors[mut_type]))
    
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
        scale_y_continuous(breaks=seq(0, ymax, scale_y)) +
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
        scale_y_continuous(breaks=seq(0, ymax, scale_y)) +
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

##
## Deprecated variants
##

#'
#' This function has been removed.  Use 'plot_profiles' instead.
#'
#' @param vcf        A GRanges object
#' @param ref_genome The reference genome
#'
#' @return Character vector with the context of the base substitutions
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' mut_context <- mut_context(vcfs[[1]], ref_genome)
#'
#' @seealso
#' \code{\link{mut_context}}
#'
#' @export

plot_96_profile <- function(vcf, ref_genome)
{
  .Defunct("plot_profiles", package="MutationalPatterns",
           msg=paste("This function has been renamed. Use",
                     "'plot_profiles' instead."))
}

