# Global package variables

# Default colours for mutation spectrum plotting
COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")

COLORS7 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#E98C7B", "#D4D2D2", "#ADCC54",
    "#F0D0CE")

COLORS10 = c(
    "#A7CBDD", "#2776AE", "#B1DE8B",
    "#389C2F", "#F49C9B", "#D41F1E",
    "#FABC73", "#FB8101", "#C3AFD5",
    "#7D599E"
)

SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)

C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one 
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
TRIPLETS_192_trans = paste(rep(TRIPLETS_96, each=2), c("transcribed", "untranscribed"), sep="-")
TRIPLETS_192_rep = paste(rep(TRIPLETS_96, each=2), c("left", "right"), sep="-")

DBS = c(
  'AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 
  'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 
  'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT',
  'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 
  'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG',
  'GC>AA', 'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 
  'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG', 'TA>GT',
  'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 
  'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT',
  'TT>AA', 'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG')

DBS_strand = rep(DBS, each = 2)
DBS_trans = paste(rep(DBS, each=2), c("transcribed", "untranscribed"), sep="-")
DBS_rep = paste(rep(DBS, each=2), c("left", "right"), sep="-")

SUBSTITUTIONS_DBS = c("AC>NN", "AT>NN", "CC>NN","CG>NN","CT>NN","GC>NN","TA>NN","TC>NN","TG>NN","TT>NN")
ALT_DBS = do.call(rbind, strsplit(DBS, ">"))[,2]

STRAND = rep(c("U","T"), 96)
STRAND_DBS = rep(c("U","T"), 78)
DNA_BASES = c("A", "C", "G", "T")
