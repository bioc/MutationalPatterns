#Script to convert a signature file from SIGNAL into the correct format for MutationalPatterns
library(dplyr)
library(stringr)

format_signatures = function(fname){
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

format_signatures("~/Downloads/snv_SIGNAL_tissuespec.txt")
format_signatures("~/Downloads/snv_SIGNAL_ref.txt")
format_signatures("~/Downloads/snv_SIGNAL_exposure.txt")

#DBS data was not on signature website, but has been extracted from the paper:
# "A Compendium of Mutational Signatures of Environmental Agents"