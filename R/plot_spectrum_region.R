#' Plot point mutation spectrum per genomic region
#' 
#' A spectrum similar to the one from 'plot_spectrum()' is plotted.
#' However the spectrum is plotted separately per genomic region.
#' As input it takes a 'type_occurrences' matrix that was calculated per genomic region.
#' To get a 'type_occurrences' matrix per region,
#' first use the 'split_muts_region()' function on a GR or GRangesList.
#' Then use the 'mut_type_occurrences' as you would normally.
#' The by, colors and legend argument work the same as in'plot_spectrum()'.
#'    
#' @param type_occurrences Type occurrences matrix
#' @param by Optional grouping variable
#' @param mode 'absolute' or 'relative'.
#' When relative, the number of variants will be shown
#' divided by the total number of variants in that sample and genomic region.
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = TRUE
#' @return Spectrum plot by genomic region
#'
#' @import ggplot2
#' @importFrom magrittr %>% 
#'
#' @examples
#' ## See the 'split_muts_region()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/grl_split_region.rds",
#'                 package="MutationalPatterns"))
#' 
#'## Or plot spectrum per tissue
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#'  ## Load a reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#'
#' ## Get the type occurrences for all VCF objects.
#' type_occurrences = mut_type_occurrences(grl, ref_genome)
#' 
#' ## Plot the point mutation spectrum per sample type and per genomic region
#' plot_spectrum_region(type_occurrences, by = tissue)
#' 
#' ## Plot the absolute point mutation spectrum per sample type and per genomic region
#' plot_spectrum_region(type_occurrences, mode = "absolute")
#' 
#' sample_names <- c(
#' "colon1", "colon2", "colon3",
#' "intestine1", "intestine2", "intestine3",
#' "liver1", "liver2", "liver3")
#' 
#' ## Plot the point mutation spectrum per individual sample and per genomic region
#' plot_spectrum_region(type_occurrences, by = sample_names)
#' 
#' 
#' @seealso
#' \code{\link{read_vcfs_as_granges}},
#' \code{\link{mut_type_occurrences}},
#' \code{\link{plot_spectrum}}
#'
#' @export
#' 
plot_spectrum_region = function(type_occurrences, by = NA, mode = c("relative", "absolute"), colors = NULL, legend = TRUE){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    `C>T at CpG` = `C>T other` = type = amount = stdev = tot_muts = lower = upper = NULL
    
    
    if (is.null(colors)) {
        colors = COLORS6
    }
    mode = match.arg(mode)
    
    row_names = rownames(type_occurrences)
    max_dots_in_name = row_names %>% 
        stringr::str_count("\\.") %>% 
        max()
    if (max_dots_in_name > 1){
        stop("The row names of the type_occurrences dataframe too many dots. 
             There should only be a dot in between the sample name and the type")
    }
    
    #Get sample names and features
    sample_names = stringr::str_remove(row_names, "\\..*")
    feature = stringr::str_remove(row_names, ".*\\.")
    feature = factor(feature, levels = unique(feature))
    
    #Remove CpG split
    type_occurrences = type_occurrences %>% 
        dplyr::select(-`C>T at CpG`, -`C>T other`)
    
    #Count total muts
    tot_muts_tb = tibble::tibble("sample" = sample_names, "tot_muts" = rowSums(type_occurrences))
    
    #Normalize if that mode is chosen
    if (mode == "relative"){
        type_occurrences = type_occurrences %>% 
            as.matrix() %>% 
            prop.table(1)
        y_lab = "Relative contribution"
    } else{
        y_lab = "Contribution"
    }
    
    #Create long tb
    tb = type_occurrences %>% 
        as.data.frame() %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(sample = sample_names, feature = feature) %>% 
        tidyr::gather(key = "type", value = "amount", -sample, -feature)
    
    
    #Combine samples
    if (is_na(by)){
        by = "all"
    }
    tb_by = tibble::tibble("sample" = unique(tb$sample),
                           "by" = by)
    tb = tb %>% 
        dplyr::left_join(tb_by, by = "sample") %>% 
        dplyr::group_by(by, feature, type) %>% 
        dplyr::summarise(stdev = sd(amount), amount = mean(amount)) %>% 
        dplyr::ungroup() %>% 
        dplyr::rename(sample = by) %>% 
        dplyr::mutate(lower = amount - stdev, upper = amount + stdev)
    
    #Create mut counts of samples
    tot_muts_tb = tot_muts_tb %>% 
        dplyr::left_join(tb_by, by = "sample") %>% 
        dplyr::group_by(by) %>% 
        dplyr::summarise(tot_muts = sum(tot_muts))
    
    #Create facet labels
    facet_labs_y = stringr::str_c(tot_muts_tb$by, " (n = ", tot_muts_tb$tot_muts, ")")
    names(facet_labs_y) = tot_muts_tb$by
    
    #Create figure
    #Suppress warning about using alpha.
    withCallingHandlers({
        fig = ggplot(tb, aes(x = type, y = amount, fill = type, alpha = feature)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black", cex = 0.5) + 
            facet_grid(. ~ sample, labeller = labeller(sample = facet_labs_y)) + 
            scale_fill_manual(values = colors) + 
            scale_alpha_discrete(range = c(1, 0.4)) + 
            labs(y = y_lab, x = "") +
            scale_x_discrete(breaks = NULL) +
            theme_bw()
    }, warning = function(w) {
        if (grepl("Using alpha for a discrete variable is not advised.", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    
    #Add errorbars
    if (sum(is.na(tb$stdev)) == 0){
        fig = fig + 
            geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge")
    }
    
    #Remove legend if required
    if (legend == FALSE)
        fig = fig + theme(legend.position="none")
    
    return(fig)
}
