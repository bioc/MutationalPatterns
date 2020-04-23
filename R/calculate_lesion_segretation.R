#' Calculate the amount of lesion segregation for a GRangesList or GRanges object.
#' 
#' The amount of lesion segregation is calculated per GRanges object.
#' The results are then combined in a table.
#'
#' @param gr GRanges object
#' @param sample_name The name of the sample
#' @param split_by_type Boolean describing whether the lesion 
#' segregation should be calculated for all SNVs together or per 96 substitution context
#' @param ref_genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome.
#'               Only needed when split_by_type is TRUE
#'
#' @return A tibble containing the amount of lesions segregation per sample
#' @importFrom magrittr %>% 
#'
#'@examples
calculate_lesion_segregation = function(grl, sample_names, split_by_type = FALSE , ref_genome = NA){
    
    if (length(grl) != length(sample_names)){
        stop("The grl and the sample_names should be equally long.", call. = F)
    }
    
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        strand_tb = purrr::map2(gr_l, sample_names, function(gr, sample_name){
            calculate_lesion_segregation_gr(gr, sample_name, split_by_type, ref_genome)
            }) %>% 
            do.call(rbind, .)
    } else if (inherits(grl, "GRanges")){
        strand_tb = calculate_lesion_segregation_gr(grl, sample_names, split_by_type, ref_genome)
    } else{
        not_gr_or_grl(grl)
    }
    strand_tb = strand_tb %>% 
        dplyr::mutate(sample_name = sample_names,
                      p_adjusted = p.adjust(p.value, method = "fdr"),
                      switching = estimate <= 0.4 & p_adjusted <= 0.05)
    return(strand_tb)
}

#' Calculate the amount of lesion segregation for a singe GRanges object.
#'
#' @param gr GRanges object
#' @param sample_name The name of the sample
#' @param split_by_type Boolean describing whether the lesion 
#' segregation should be calculated for all SNVs together or per 96 substitution context.
#' @param ref_genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome.
#'               Only needed when split_by_type is TRUE
#' @return A tibble containing the amount of lesions segregation for a single sample
#' @importFrom magrittr %>% 
#'
calculate_lesion_segregation_gr = function(gr, sample_name = "sample", split_by_type = FALSE, ref_genome = NA){
    
    
    if (!length(gr)){
        message(str_c("No mutations present in sample: ", sample_name,
                      "\n Returning NA"))
        return(NA)
    }
    
    gr = get_strandedness_gr(gr)
    tb = get_strandedness_tb(gr)
    
    if (split_by_type){
        
        #Split gr according to the 96 substitution context.
        cnd = tryCatch(suppressWarnings({seqlevelsStyle(gr) = "UCSC"}),
                       error = function(cnd) cnd)
        if (inherits(cnd, "error")){
            message(str_c("Could not change seqlevelstyle in sample: ", sample_name, ".",
                          "\n Returning NA"))
            return(NA)
        }
        type_context = type_context(gr, get(ref_genome))
        full_context = stringr::str_c(substr(type_context$context, 1, 1), 
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
        message(str_c("No multiple mutations in one chromosome with context present in sample: ", sample_name,
                      "\n Returning NA"))
        return(NA)
    }
    
    #Calculate if the number of strand switches is significantly different from expected.
    stat_tb = binom.test(x = res$x, n = res$n, p = 0.5) %>% 
        broom::tidy()
    return(stat_tb)
}

#' Determine the strands of a GRanges object
#'
#' @param gr A GRanges object
#'
#' @return A GRanges object where the strands have been set.
#' 
get_strandedness_gr = function(gr){
    check_no_indels(gr)
    strand(gr) = ifelse(as.vector(gr$REF) %in% c("C", "T"), "+", "-")
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
calculate_strand_switches = function(tb){
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
#' @importFrom magrittr %>%  
calculate_strand_switch = function(strand){
    switches = strand != dplyr::lead(strand)
    switches = switches %>% 
        stats::na.omit() %>% 
        as.vector()
    return(switches)
}