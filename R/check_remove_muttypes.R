#' Checks that there are no variants with multiple alt alleles.
#'
#' @details 
#' This function checks for a GRanges/GRangesList object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_multi_alts = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr = unlist(grl)
    } else if (inherits(grl, "GRanges")){
        gr = grl
    } else{
        not_gr_or_grl(grl)
    }
    check_no_multi_alts_gr(gr)
    invisible(grl)
}

#' Checks that there are no variants with multiple alt alleles.
#'
#' @details 
#' This function checks for a single GRanges object, that there are no variants with multiple alternative alleles.
#' It works by checking that the number of variants is equal to the number of alternative alleles.
#' It trows an error when this is not the case.
#' 
#' @param gr GRanges object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_multi_alts_gr = function(gr){
    alt = gr$ALT
    nr_alts = alt %>% 
        unlist() %>% 
        length()
    if (length(gr) != nr_alts){
        stop("There should not be any variants with multiple alternative alleles. Please remove these variants with remove_multi_alts_variants.", call. = F)
    }
    invisible(gr)
}

#' Checks that there are no SNVs/MNVs.
#'
#' @details 
#' This function checks for a GRanges/GRangesList object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#' 
#' @param grl GRanges/GrangesList object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 

check_no_snvs = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr = unlist(grl)
    } else if (inherits(grl, "GRanges")){
        gr = grl
    } else{
        not_gr_or_grl(grl)
    }
    check_no_snvs_gr(gr)
    invisible(grl)
}

#' Checks that there are no SNVs/MNVs.
#'
#' @details 
#' This function checks for a single GRanges object, that there are no SNVs/MNVs.
#' It trows an error when this is not the case.
#' 
#' @param gr GRanges object
#' 
#' @return Invisibly returns its input when it succeeds.
#' @examples
#' 
#' @importFrom magrittr %>%

check_no_snvs_gr = function(gr){
    snv_f = find_snv(gr)
    nr_snv = sum(snv_f)
    snv_present = nr_snv >= 1
    if (snv_present){
        stop(stringr::str_c("There seem to be ", nr_snv, " SNVs present in your data. Make sure to remove all SNVs with the remove_snvs function."), call. = F)
    }
    invisible(gr)
}

#' Removes variants with multiple alt alleles.
#'
#' @details 
#' This function removes variants with multiple alternative alleles for a GRanges/GRangesList object.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return A filtered version of the input GRanges/GRangesList object.
#' @examples
#' 
#' @export

remove_multi_alts_variants = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        grl = purrr::map(gr_l, remove_multi_alts_variants_gr) %>% 
            GRangesList()
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = remove_multi_alts_variants_gr(grl)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}

#' Removes variants with multiple alt alleles.
#'
#' @details 
#' This function removes variants with multiple alternative alleles for a single GRanges object.
#' 
#' @param gr GRanges object
#' 
#' @return A filtered version of the input GRanges object.
#' @examples
#' 

remove_multi_alts_variants_gr = function(gr){
    alt = gr$ALT
    gr = gr[elementNROWS(alt) == 1]
    return(gr)
}

#' Removes SNV/MNV variants.
#'
#' @details 
#' This function removes SNV/MNV variants for a GRanges/GrangesList object.
#' 
#' @param grl GRanges/GRangesList object
#' 
#' @return A filtered version of the input GRanges/GrangesList object.
#' @examples
#' 
#' @export

remove_snvs = function(grl){
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        grl = purrr::map(gr_l, remove_snvs_gr) %>% 
            GRangesList()
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = remove_snvs_gr(grl)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}

#' Removes SNV/MNV variants.
#'
#' @details 
#' This function removes SNV/MNV variants for a single GRanges object.
#' 
#' @param gr GRanges object
#' 
#' @return A filtered version of the input GRanges object.
#' @examples
#' 

remove_snvs_gr = function(gr){
    snv_f = find_snv(gr)
    gr = gr[!snv_f]
    return(gr)
}

#' Identifies SNVs/MNVs
#'
#' @details 
#' This function finds SNVs/MNVs for a single GRanges object.
#' 
#' 
#' @param gr GRanges object
#' 
#' @return A boolean vector. It's TRUE for SNVs/MNVs and FALSE for Indels
#' @examples
#' 

find_snv = function(gr){
    check_no_multi_alts(gr)
    snv_f = width(gr$REF) == width(unlist(gr$ALT))
    return(snv_f)
}

#' Throw an error for when an argument should be a GRanges/GRangesList object.
#'
#' @details 
#' This function is called by other functions, when an argument should be a GRanges/GrangesList object, but isn't.
#' It throws an error showing the actual class of the argument.
#' 
#' @param arg Argument. Any object that should have been a GRanges/GrangesList, but isn't.
#' 
#' @examples
#' \dontrun{
#' a = 1
#' not_gr_or_grl(a)
#' }
#' 

not_gr_or_grl = function(arg){
    arg_name = deparse(substitute(arg))
    arg_class = class(arg)[[1]]
    stop(stringr::str_c(arg_name, " should be a CompressedGRangesList or GRanges object, instead it is a ", arg_class, " object.
               Please provide an object of the correct class."))
}

#' Check that the seqnames of a GRanges object are present in a ref_genome
#'
#' @details 
#' This function tests that the variants in a GRanges object have seqnames, that match the supplied ref_genome.
#' If this is not the case an error is thrown, with suggestions how to fix this.
#' 
#' @param gr GRanges object
#' @param ref_genome BSGenome reference genome object
#' 
#' @examples
#' 

check_chroms = function(gr, ref_genome){
    gr_seqnames = as.vector(seqnames(gr))
    ref_seqnames = seqnames(get(ref_genome))
    shared_chroms = intersect(gr_seqnames, ref_seqnames)
    if (!length(shared_chroms)){
        stop("The input GRanges and the ref_genome share no seqnames (chromosome names). 
             Do they use the same seqlevelsStyle? An example of how to fix this is show below.
             You can change the seqlevelStyle with: `seqlevelsStyle(grl) = 'UCSC'", call. = F)
    }
    gr_seqlevels = levels(seqnames(gr))
    gr_seqlevels_notref = unique(gr_seqlevels[!gr_seqlevels %in% ref_seqnames])
    gr_seqlevels_notref = stringr::str_c(gr_seqlevels_notref, collapse = ", ")
    if(length(gr_seqlevels_notref)){
        stop(stringr::str_c("The following seqlevels (chromosome names) occur in the input GRanges,
        but are not present in the ref_genome: ", gr_seqlevels_notref, ". An example of how to fix this is show below.
        First select the chromosomes you want to keep with: 
        `chromosomes = paste0('chr', c(1:22,'X'))`
        You can then remove variants in other chromosomes with: 
        `seqlevels(grl, pruning.mode = 'tidy') = chromosomes`"), call. = F)
    }
    invisible(gr)
}
