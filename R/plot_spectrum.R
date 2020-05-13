#' Plot point mutation spectrum
#'    
#' @param type_occurrences Type occurrences matrix
#' @param CT Distinction between C>T at CpG and C>T at other
#' sites, default = FALSE
#' @param by Optional grouping variable
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = TRUE
#' @return Spectrum plot
#'
#' @import ggplot2
#' @importFrom magrittr %>% 
#'
#' @examples
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#' 
#'
#' ## Load a reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences = mut_type_occurrences(vcfs, ref_genome)
#' 
#' ## Plot the point mutation spectrum over all samples
#' plot_spectrum(type_occurrences)
#'
#' ## Or with distinction of C>T at CpG sites
#' plot_spectrum(type_occurrences, CT = TRUE)
#'
#' ## Or without legend
#' plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
#'
#' ## Or plot spectrum per tissue
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' plot_spectrum(type_occurrences, by = tissue, CT = TRUE)
#'
#' ## You can also set custom colors.
#' my_colors = c("pink", "orange", "blue", "lightblue",
#'                 "green", "red", "purple")
#'
#' ## And use them in a plot.
#' plot_spectrum(type_occurrences,
#'                 CT = TRUE,
#'                 legend = TRUE,
#'                 colors = my_colors)
#'
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_type_occurrences}}
#'
#' @export

plot_spectrum = function(type_occurrences, CT=FALSE, by = NA, colors = NA, legend = TRUE)
{
  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
    value = nmuts = sub_type = variable = error_pos = stdev = NULL


    # If colors parameter not provided, set to default colors
    if (is_na(colors))
        colors = COLORS7

    # Check color vector length
    if (length(colors) != 7)
        stop("Colors parameter: supply color vector with length 7")

    # Distinction between C>T at CpG or not
    if (CT == FALSE)
        type_occurrences = type_occurrences[,1:6] 
    else
        type_occurrences = type_occurrences[,c(1:2,8,7,4:6)]

    # If grouping variable not provided, set to "all"
    if (is_na(by))
        by="all"

    # Reshape the type_occurences for the plotting
    tb = type_occurrences %>% 
      tibble::rownames_to_column("sample") %>% 
      dplyr::mutate(by = by) %>% #Add user defined grouping
      tidyr::pivot_longer(c(-sample, -by), names_to = "variable", values_to = "nmuts") %>% #Make long format
      dplyr::group_by(sample) %>% 
      dplyr::mutate(value = nmuts / sum(nmuts)) %>% #Calculate relative values
      dplyr::ungroup() %>% 
      dplyr::group_by(by, variable) %>% #Summarise per group and mutation type
      dplyr::summarise(mean = mean(value), stdev = stats::sd(value),
                       total_individuals = sum(value), total_mutations = sum(nmuts)) %>% 
      dplyr::mutate(total_individuals = sum(total_individuals), total_mutations = sum(total_mutations)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(total_mutations = prettyNum(total_mutations, big.mark = ","), #Make pretty and add subtypes
                    total_mutations = paste("No. mutations = ", total_mutations),
                    sub_type = stringr::str_remove(variable, " .*"),
                    variable = factor(variable, levels = unique(variable)),
                    error_pos = mean)
    
    # Define colors for plotting
    if(CT == FALSE)
        colors = colors[c(1,2,3,5:7)]

    # C>T stacked bar (distinction between CpG sites and other)
    else
    {
        # Adjust positioning of error bars for stacked bars
        # mean of C>T at CpG should be plus the mean of C>T at other
        CpG = which(tb$variable == "C>T at CpG")
        other = which(tb$variable == "C>T other")
        tb$error_pos[CpG] = tb$error_pos[other] + tb$error_pos[CpG]
    }

    # Make barplot
    plot = ggplot(data=tb, aes(x=sub_type,
                                y=mean,
                                fill=variable,
                                group=sub_type)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values=colors, name="Point mutation type") +
        theme_bw() +
        xlab("") +
        ylab("Relative contribution") +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                panel.grid.major.x = element_blank())
    
    # check if standard deviation error bars can be plotted
    if(sum(is.na(x$stdev)) > 0)
      warning("No standard deviation error bars can be plotted, because there is only one sample per mutation spectrum")
    else
      plot = plot + geom_errorbar(aes(ymin=error_pos-stdev, 
                                      ymax=error_pos+stdev), width=0.2)
    
      
    # Facetting
    if (length(by) == 1)
        plot = plot + facet_wrap( ~ total_mutations)
    else
        plot = plot + facet_wrap(by ~ total_mutations)

    # Legend
    if (legend == FALSE)
        plot = plot + theme(legend.position="none")

    return(plot)
} 
