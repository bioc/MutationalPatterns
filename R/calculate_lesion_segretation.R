#' Calculate the amount of lesion segregation for a GRangesList or GRanges object.
#' 
#' The amount of lesion segregation is calculated per GRanges object.
#' The results are then combined in a table.
#' It's possible to calculate the lesion segregation separately per 96 substitution context,
#' when using the binomial test. The results are then automatically added back up together.
#'
#' @param grl GRangesList or GRanges object
#' @param sample_names The name of the sample
#' @param test The statistical test that should be used. Possible values:
#'              * 'binomial' Binomial test based on the number of strand switches. (Default);
#'              * 'walf-wolfowitz' Statistical test that checks if the strands are randomly distibuted.;
#' @param split_by_type Boolean describing whether the lesion 
#' segregation should be calculated for all SNVs together or per 96 substitution context. (Default: FALSE)
#' @param ref_genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome.
#'               Only needed when split_by_type is TRUE
#'
#' @return A tibble containing the amount of lesions segregation per sample
#' @importFrom magrittr %>% 
#' @seealso
#' \code{\link{plot_lesion_segregation}}
#' @family Lesion_segregation
#' @export
#' @examples
#'
#' ## See the 'read_vcfs_as_granges()' example for how we obtained the
#' ## following data:
#' grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
#'                 package="MutationalPatterns"))
#' 
#' ## Set the sample names
#' sample_names <- c(
#'    "colon1", "colon2", "colon3",
#'    "intestine1", "intestine2", "intestine3",
#'    "liver1", "liver2", "liver3")
#'                                  
#' ## Load the corresponding reference genome.
#' ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#' 
#' ## Calculate lesion segregation
#' lesion_segretation = calculate_lesion_segregation(grl, sample_names)
#' 
#' ## Calculate lesion segregation per 96 base type
#' lesion_segretation_by_type = calculate_lesion_segregation(grl, sample_names, 
#' split_by_type = TRUE, ref_genome = ref_genome)
#' 
#' ## Calculate lesion segregation using the walf-wolfowitz test.
#' lesion_segregation_walf = calculate_lesion_segregation(grl, 
#'                                                        sample_names, 
#'                                                        test = "walf-wolfowitz")
#' 
calculate_lesion_segregation = function(grl, 
                                        sample_names, 
                                        test = c("binomial", "walf-wolfowitz"),
                                        split_by_type = FALSE, 
                                        ref_genome = NA){
    
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    p.value = . = sample_name = NULL
    
    #Validate arguments
    test = match.arg(test)
    if (test == "walf-wolfowitz" & split_by_type){
        stop("The 'split_by_type' argument can only be used with the binomial test",
             call. = F)
    }
    if (length(grl) != length(sample_names)){
        stop("The grl and the sample_names should be equally long.", call. = F)
    }
    
    if (split_by_type){
        if (is_na(ref_genome)){
            stop("The ref_genome needs to be set when split_by_type = TRUE")
        }
    }
    
    #Perform lesion segregation on each GR
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        strand_tb = purrr::map2(gr_l, sample_names, function(gr, sample_name){
            calculate_lesion_segregation_gr(gr, sample_name, test, split_by_type, ref_genome)
            }) %>% 
            do.call(rbind, .)
    } else if (inherits(grl, "GRanges")){
        strand_tb = calculate_lesion_segregation_gr(grl, sample_names, test, split_by_type, ref_genome)
    } else{
        not_gr_or_grl(grl)
    }
    
    #Add final columns to output
    strand_tb = strand_tb %>% 
        dplyr::mutate(sample_name = sample_names,
                      fdr = p.adjust(p.value, method = "fdr")) %>% 
        dplyr::select(sample_name, dplyr::everything())
    return(strand_tb)
}

#' Calculate the amount of lesion segregation for a singe GRanges object.
#'
#' @param gr GRanges object
#' @param sample_name The name of the sample
#' @param test The statistical test that should be used. Possible values:
#'              * 'binomial' Binomial test based on the number of strand switches. (Default);
#'              * 'walf-wolfowitz' Statistical test that checks if the strands are randomly distibuted.;
#' @param split_by_type Boolean describing whether the lesion 
#' segregation should be calculated for all SNVs together or per 96 substitution context.
#' @param ref_genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome.
#'               Only needed when split_by_type is TRUE
#' @return A tibble containing the amount of lesions segregation for a single sample
#' @noRd
#'
calculate_lesion_segregation_gr = function(gr, 
                                           sample_name = "sample", 
                                           test = c("binomial", "walf-wolfowitz"),
                                           split_by_type = FALSE, 
                                           ref_genome = NA){
    
    
    #Check if mutations are present.
    if (!length(gr)){
        message(paste0("No mutations present in sample: ", sample_name,
                      "\n Returning NA"))
        return(NA)
    }
    
    #Get strand info
    gr = get_strandedness_gr(gr)
    tb = get_strandedness_tb(gr)
    
    if (test == "binomial"){
        #Perform analysis per base substitution type
        if (split_by_type){
            
            
            #Split gr according to the 96 substitution context.
            cnd = tryCatch(suppressWarnings({GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"}),
                           error = function(cnd) cnd)
            if (inherits(cnd, "error")){
                message(paste0("Could not change seqlevelstyle in sample: ", sample_name, ".",
                              "\n Returning NA"))
                return(NA)
            }
            check_chroms(gr, ref_genome)
            type_context = type_context(gr, ref_genome)
            full_context = paste0(substr(type_context$context, 1, 1), 
                                 "[", type_context$types, "]", 
                                 substr(type_context$context, 3, 3))
            tb_l = split(tb, full_context)
            
            #Calculate strand switches for each of the 96 substitutions.
            res_l = purrr::map(tb_l, calculate_strand_switches)
            x = purrr::map(res_l, "x") %>% 
                unlist() %>% 
                sum()
            n = purrr::map(res_l, "n") %>% 
                unlist() %>% 
                sum()
            res = list("x" = x, "n" = n)
        } else{
            #Calculate strand switches
            res = calculate_strand_switches(tb)
        }
        
        #Check if mutations are present
        if (res$n == 0){
            message(paste0("No multiple mutations in one chromosome with context present in sample: ", sample_name,
                          "\n Returning NA"))
            return(NA)
        }
        
        #Calculate if the number of strand switches is significantly different from expected.
        res = binom.test(x = res$x, n = res$n, p = 0.5)
        
        #Add all results together in a tibble
        stat_tb = tibble::tibble(p.value = res$p.value, 
                                 percent_strand_switches = res$estimate,
                                 conf_low = res$conf.int[[1]], 
                                 conf_high = res$conf.int[[2]], 
                                 nr_strand_switches = res$statistic, 
                                 max_possible_switches = res$parameter)
    } else{
        #calculate if there is a significant deviation using the walf_wolfowitz_test
        wolfowitz = walf_wolfowitz_test(tb$strand)
        stat_tb = tibble::tibble(p.value = wolfowitz$p,
                                 sd = wolfowitz$sd,
                                 nr_total_runs = wolfowitz$runs_total)
    }
    
    return(stat_tb)
}

#' Determine the strands of a GRanges object
#'
#' @param gr A GRanges object
#'
#' @return A GRanges object where the strands have been set.
#' @noRd
#' 
get_strandedness_gr = function(gr){
    check_no_indels(gr)
    strand(gr) = ifelse(as.vector(get_ref(gr)) %in% c("C", "T"), "+", "-")
    GenomeInfoDb::seqlevelsStyle(gr) = "NCBI" #This takes less space when plotting
    
    if (length(gr)){
        GenomeInfoDb::seqlevels(gr) = GenomeInfoDb::seqlevelsInUse(gr)
    }
    return(gr)
}

#' Convert a GRanges object with strand info to a tibble
#'
#' @param gr A GRanges object where the strands have been set.
#'
#' @return A tibble with strand information
#' @importFrom magrittr %>% 
#' @noRd
#' 
get_strandedness_tb = function(gr){
    tb = as.data.frame(gr) %>%
        tibble::as_tibble() %>% 
        dplyr::mutate(strand = droplevels(strand),
                      y = dplyr::recode(strand, "+" = 1, "-" = 0),
                      start_mb = start/1000000)
    return(tb)
}

#' Calculate the total number of variants and strand switches.
#'
#' @param tb A tibble with strand information
#'
#' @return A list containing the total number of variants and the number of strand switches
#' @importFrom magrittr %>%  
#' @noRd
#' 
calculate_strand_switches = function(tb){
    # These variables use non standard evaluation.
    # To avoid R CMD check complaints we initialize them to NULL.
    . = NULL
    
    strands_l = split(tb$strand, tb$seqnames)
    switches = purrr::map(strands_l, calculate_strand_switch) %>% 
        do.call("c", .)
    res = list("x" = sum(switches), "n" = length(switches))
    return(res)
}

#' Calculate the number of strand switches in a single chromosome.
#'
#' @param strand A vector containing the strand of variants
#'
#' @return A boolean vector describing for each variant 
#' if it switched strands with the previous variant.
#' @noRd
#' @importFrom magrittr %>%  
#' 
calculate_strand_switch = function(strand){
    switches = strand != dplyr::lead(strand)
    switches = switches %>% 
        stats::na.omit() %>% 
        as.vector()
    return(switches)
}


#' Perform the walf_wofowitz test for strands.
#' 
#' This statistical test, tests whether each element in the sequence is 
#' independently drawn from the same distribution.
#'
#' @param strands A vector of strands
#'
#' @return a p value
#'
#' @noRd
#' 
walf_wolfowitz_test = function(strands){
    
    #Remove factor
    strands = as.character(strands)
    
    #Determine sizes
    n1 = sum(strands == "+")
    n2 = sum(strands == "-")
    n = n1 + n2
    
    #Determine number of + and - runs
    runs = rle(strands)
    r1 = length(runs$lengths[runs$values=="+"])
    r2 = length(runs$lengths[runs$values=="-"])  
    
    #Calculate total number of runs
    runs_total = r1+r2
    
    #Calculate mean
    mean_val = 2*n1*n2/(n) + 1
    
    #Calculate variance and sd
    variance = (mean_val-1)*(mean_val-2)/(n-1)
    sd = sqrt(variance)
    
    #Calculate p value
    p <- stats::pnorm((runs_total - mean_val) / sd)
    
    #Make two-sided
    p <- 2*min(p,1-p)
    
    return(list("p" = p, "sd" = sd, "runs_total" = runs_total))
}
