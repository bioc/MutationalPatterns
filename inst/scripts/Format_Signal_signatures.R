#Script to convert a signature file from SIGNAL into the correct format for MutationalPatterns
library(dplyr)
library(stringr)
library(readr)

format_SIGNAL_signatures = function(fname){
    signatures =read.table(fname, 
                           header = TRUE, 
                           sep = "\t", 
                           stringsAsFactors = FALSE)
    
    colnames(signatures)[1] = "Type_subtype"
    signatures = signatures %>% 
        dplyr::mutate(Type = str_replace(Type_subtype, ".*\\[(.*)\\].*", "\\1"),
                      SubType = str_remove_all(Type_subtype, "\\[|\\]|\\>[A-Z]")) %>% 
        dplyr::select(-Type_subtype) %>% 
        dplyr::select(Type, SubType, everything())
    
    fname_base = basename(fname)
    out_path = file.path("inst", "extdata", "signatures", fname_base)
    write.table(signatures, 
                out_path, 
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    invisible(0)
}

format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_tissue.txt")
format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_reference.txt")
format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_exposure.txt")

#DBS data was not on signature website, but has been extracted from the paper:
# "A Compendium of Mutational Signatures of Environmental Agents

#Add sparse signatures from the paper:
# "De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures"
signatures = read_tsv("~/Downloads/snv_SPARSE.txt", 
                      col_types = cols(.default = "d", sig = "c"), 
                      locale=locale(decimal_mark = ","))
signatures = signatures %>% 
    tidyr::pivot_longer(-sig, names_to = "Type_subtype", values_to = "values") %>% 
    tidyr::pivot_wider(names_from = sig, values_from = values) %>% 
    dplyr::mutate(Type = str_replace(Type_subtype, ".*\\[(.*)\\].*", "\\1"),
                  SubType = str_remove_all(Type_subtype, "\\[|\\]|\\>[A-Z]")) %>% 
    dplyr::select(-Type_subtype) %>% 
    dplyr::select(Type, SubType, everything())

write.table(signatures, 
            "inst/extdata/signatures/snv_SPARSE_reference.txt", 
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

