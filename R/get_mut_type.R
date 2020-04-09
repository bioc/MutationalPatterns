#' Get variants with mut_type from GRanges
#' 
#' Get the variants of a certain mutation type from a GRanges or GRangesList object.
#' All other variants will be filtered out.
#' It is assumed that DBS/MBSs are called as separate SNVs.
#' They are merged into single variants.
#' The type of variant can be chosen with type.
#'
#' @param gr GRanges/GRangesList
#' @param type The type of variant that will be returned.
#'
#' @return GRanges/GrangesList of the desired mutation type.
#' 
#' @examples 
#' ## Get a GRanges list object.
#' ## See 'read_vcfs_as_granges' for more info how to do this.
#' grl <- readRDS(system.file("states/blood_grl.rds",
#'                 package="MutationalPatterns"))
#' 
#' ## Get a specific mutation type.
#' grl_snv = get_mut_type(grl, "snv")
#' grl_indel = get_mut_type(grl, "indel")
#' grl_dbs = get_mut_type(grl, "dbs")
#' grl_mbs = get_mut_type(grl, "mbs")
#' 
#' @importFrom magrittr %>%
#' @export
get_mut_type = function(grl, type = c("snv", "indel", "dbs", "mbs")){
    type = match.arg(type)
    
    if (inherits(grl, "CompressedGRangesList")){
        gr_l = as.list(grl)
        grl = purrr::map(gr_l, get_mut_type_gr, type) %>% 
            GenomicRanges::GRangesList()
        return(grl)
    } else if (inherits(grl, "GRanges")){
        gr = get_mut_type_gr(grl, type)
        return(gr)
    } else{
        not_gr_or_grl(grl)
    }
}


#' Get variants with mut_type from GRanges
#' 
#' Get the variants of a certain mutation type from a GRanges object.
#' All other variants will be filtered out.
#' It is assumed that DBS/MBSs are called as separate SNVs.
#' They are merged into single variants.
#' The type of variant can be chosen with type.
#'
#' @param gr GRanges
#' @param type The type of variant that will be returned.
#'
#' @return GRanges of the desired mutation type.

get_mut_type_gr = function(gr, type = c("snv", "indel", "dbs", "mbs")){
    type = match.arg(type)
    
    #Filter out bad variants
    gr = remove_multi_alts_variants(gr)
    
    #Indels
    if (type == "indel"){
        gr = remove_snvs(gr)
        return(gr)
    }
    
    #Split SNVs into SNVs/DBS/MNVs
    gr = remove_indels(gr)
    gr_l = split_mbs_gr(gr)
    
    if (type == "snv"){
        gr = gr_l$`1`
    } else if (type == "dbs"){
        gr = gr_l$`2`
    } else if (type == "mbs"){
        not_mbs_f = names(gr_l) %in% c(1,2) #Determine which elements of the list are MNVs
        if (sum(!not_mbs_f) == 0){
            gr = GenomicRanges::GRanges()
        } else{
            gr_l = gr_l[!not_mbs_f]
            gr = unlist(GenomicRanges::GRangesList(gr_l))
        }
        if (is.null(gr)){
            gr = GenomicRanges::GRanges()
        }
    }
    return(gr)
}

#' Split SNV/MNV into SNV, DBS and MNVs of different sizes
#'
#'A function that splits a GRanges object into a list of GRanges.
#'Mutations are split based on whether they are a SNV/MNV. They are split based on their size.
#'
#'
#'This function assumes that a DBS/MNV is called as separate SNVs, that are located next to each other.
#'This function works by checking the location of these SNVs. 
#'If multiple SNVs are located next to eachother they are considered a DBS/MNV.
#'A CAT>TAG would be considered as a MNV of length 3.
#'The largest variants are then selected and put into a separate GRanges object.
#'This is done recursively untill all variants are separated based on their size.
#'This way the CAT>TAG mutation will not be seen as two separate DBS variants.
#'This results in a list of granges.
#'If merge_muts is TRUE, then the SNVs that make up a DBS/MNV will be merged.
#'
#'
#' @param gr A GRanges object
#' @param merge_muts Boolean value. If TRUE, the mutations are merged
#'
#' @return A list of granges
#' @importFrom magrittr %>% 
split_mbs_gr = function(gr, merge_muts = T){
    
    check_no_indels(gr)
    
    gr = BiocGenerics::sort(gr)
    
    #Identify location of each mut and its subsequent mut.
    chroms = GenomeInfoDb::seqnames(gr) %>% 
        as.vector()
    first_chrom = chroms[-length(chroms)]
    second_chrom = chroms[-1]
    
    coords = BiocGenerics::start(gr)
    first_coords = coords[-length(coords)]
    second_coords = coords[-1]
    
    #Compare location of mut and subsequent mut. If they are only one base apart then they are sequential.
    #Sequential forms a boolean vector. A 1 means the mutations is sequential to the previous one
    sequential = first_chrom == second_chrom & first_coords + 1 == second_coords
    sequential = c(0, sequential)
    
    #Sequence is a numeric vector. A n means that a mutations is the Nth sequential base in a mutations. 
    #Example: 0100123. Here the second mutation is sequential to the first. So the first two muts are a dbs.
    #Next is a sbs. Then there is another 0 and then 3 sequential muts. This means that together they form a 4 base substitution.
    seq_rle = rle(sequential)
    sequence = sequence(seq_rle$lengths)
    sequence[sequential == 0] = 0
    
    #Find the maximum number of sequential mutations. (This is one less than the final mut size. So 1 for a dbs.)
    max_seq_l = max(sequence)
    mut_l = max_seq_l + 1
    
    #Find the locations of the variants with this mut size
    end_i_muts = which(max_seq_l == sequence)
    start_i_muts = end_i_muts - max_seq_l
    full_muts_i_l = purrr::map2(start_i_muts, end_i_muts, seq)
    full_muts_i = unlist(full_muts_i_l)
    
    #Merge dbs and mbs if the used wants this.
    if (mut_l != 1 & merge_muts == T){
        gr_sub = purrr::map(full_muts_i_l, function(i_v) merge_muts(gr[i_v])) %>% 
            do.call(base::c, .)
    } else{
        gr_sub = gr[full_muts_i]
    }
    
    #Put the mutations with the largest mutsize in a list
    gr_l = list(gr_sub)
    names(gr_l) = mut_l
    
    #Recurse on the rest mutations if there are still smaller mutions left.
    gr = gr[-full_muts_i]
    if(length(gr) == 0){
        return(gr_l)
    } else{
        gr_l_deeper = split_mbs_gr(gr)
        gr_l = c(gr_l, gr_l_deeper)
    }
    return(gr_l)
}



#' Merge variants in a GRanges object.
#' 
#' Merge all variants in a GRanges object into a single variant.
#' This is usefull when a DBS/MNV is called as separate SNVs.
#'
#' @param gr GRanges object
#'
#' @return GRanges object with a single variant.
merge_muts = function(gr){
    
    #Check input gr
    check_no_indels(gr)
    
    #Use position and other values of the first mut.
    gr_new = gr[1]
    
    #Combine refs, alts and quals
    gr_new$REF = gr$REF %>% 
        as.vector() %>% 
        stringr::str_c(collapse = "") %>%
        Biostrings::DNAStringSet()
    gr_new$ALT = gr$ALT %>%
        unlist() %>% 
        as.vector() %>% 
        stringr::str_c(collapse = "") %>%
        Biostrings::DNAStringSetList()
    gr_new$QUAL = mean(gr$QUAL, na.rm = T)
    
    #Return the new gr
    names(gr_new) = ""
    return(gr_new)
}