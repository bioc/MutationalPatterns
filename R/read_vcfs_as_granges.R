#' Read VCF files into a GRangesList
#'
#' This function reads Variant Call Format (VCF) files into a GRanges object
#' and combines them in a GRangesList.  In addition to loading the files, this
#' function applies the same seqlevel style to the GRanges objects as the
#' reference genome passed in the 'genome' parameter.
#'
#' @param vcf_files Character vector of VCF file names
#' @param sample_names Character vector of sample names
#' @param genome A string matching the name of a BSgenome library
#'               corresponding to the reference genome of your VCFs
#' @param group Selector for a seqlevel group.  All seqlevels outside
#'              of this group will be removed.  Possible values:
#'              * 'all' for all chromosomes;
#'              * 'auto' for autosomal chromosomes;
#'              * 'sex' for sex chromosomes;
#'              * 'auto+sex' for autosomal + sex chromosomes (default);
#'              * 'circular' for circular chromosomes;
#'              * 'none' for no filtering, which results in keeping all
#'                seqlevels from the VCF file.
#' @param type The mutation type that will be loaded. All other variants will be filtered out.
#'              Possible values:
#'              * 'snv'
#'              * 'indel'
#'              * 'dbs'
#'              * 'mbs'
#'              * 'all'
#' 
#' @return A GRangesList containing the GRanges obtained from 'vcf_files'
#'
#' @importFrom magrittr %>% 
#'
#' @examples
#' # The example data set consists of three colon samples, three intestine
#' # samples and three liver samples.  So, to map each file to its appropriate
#' # sample name, we create a vector containing the sample names:
#' sample_names <- c ( "colon1", "colon2", "colon3",
#'                     "intestine1", "intestine2", "intestine3",
#'                     "liver1", "liver2", "liver3" )
#'
#' # We assemble a list of files we want to load.  These files match the
#' # sample names defined above.
#' vcf_files <- list.files(system.file("extdata", 
#'                                     package="MutationalPatterns"),
#'                                     pattern = "sample.vcf", full.names = TRUE)
#'
#' # Get a reference genome BSgenome object.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library("BSgenome")
#' library(ref_genome, character.only = TRUE)
#'
#' # This function loads the files as GRanges objects.
#' # For backwards compatability reasons it only loads SNVs by default
#' vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
#' 
#' #To load all variant types use:
#' vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")
#' 
#' @export
read_vcfs_as_granges <- function(vcf_files, 
                                 sample_names, 
                                 genome, 
                                 group = "auto+sex", 
                                 type = c("snv", "indel", "dbs", "mbs", "all")){
    
    #Match argument
    type = match.arg(type)
    
    # Check sample names
    if (length(vcf_files) != length(sample_names))
        stop("Please provide the same number of sample names as VCF files")
    
    # Check the class of the reference genome
    genome <- base::get(genome)
    if (!inherits(genome, "BSgenome")){
        stop("Please provide the name of a BSgenome object.")
    }
    
    #Read vcfs
    grl <- purrr::map(vcf_files, read_single_vcf_as_grange, genome, group) %>% 
        GenomicRanges::GRangesList()
    
    #Filter for mutation type
    if (type != "all"){
        grl = get_mut_type(grl, type)
    }
    
    #Set the provided names for the samples.
    names(grl) <- sample_names
    
    return(grl)
}

#' Read a single VCF file into a GRanges object
#'
#' This function reads a Variant Call Format (VCF) file into a GRanges object
#' In addition to loading the files, this
#' function applies the same seqlevel style to the GRanges objects as the
#' reference genome passed in the 'genome' parameter.
#'
#' @param vcf_file A VCF file name
#' @param genome BSgenome object
#' @param group Selector for a seqlevel group.  All seqlevels outside
#'              of this group will be removed.  Possible values:
#'              * 'all' for all chromosomes;
#'              * 'auto' for autosomal chromosomes;
#'              * 'sex' for sex chromosomes;
#'              * 'auto+sex' for autosomal + sex chromosomes (default);
#'              * 'circular' for circular chromosomes;
#'              * 'none' for no filtering, which results in keeping all
#'                seqlevels from the VCF file.
#' @return A GRanges object
#'
read_single_vcf_as_grange = function(vcf_file, genome, group){
    
    # Use VariantAnnotation's readVcf, but only store the
    # GRanges information in memory.  This speeds up the
    # loading significantly.
    genome_name <- GenomeInfoDb::genome(genome)[[1]]
    vcf <- SummarizedExperiment::rowRanges(VariantAnnotation::readVcf(vcf_file, genome_name))
    
    # Convert to a single chromosome naming standard.
    seqlevelsStyle(vcf) <- GenomeInfoDb::seqlevelsStyle(genome)[1]
    
    #Filter for variants with the correct seqlevels
    vcf = filter_seqlevels(vcf, group, genome)
    
    return(vcf)
}

#' Filter a GRanges object based on seqlevels
#'
#' This function filters a GRanges object based on a group of seqnames.
#'
#' @param gr GRanges object 
#' @param group Selector for a seqlevel group.  All seqlevels outside
#'              of this group will be removed.  Possible values:
#'              * 'all' for all chromosomes;
#'              * 'auto' for autosomal chromosomes;
#'              * 'sex' for sex chromosomes;
#'              * 'auto+sex' for autosomal + sex chromosomes (default);
#'              * 'circular' for circular chromosomes;
#'              * 'none' for no filtering, which results in keeping all
#'                seqlevels from the VCF file.
#' @param genome BSgenome object
#' @return A GRanges object
#'
filter_seqlevels = function(gr, group, genome){
    
    groups <- c()
    if (group != "none"){
        #These variables are needed to extract the possible seqlevels
        ref_style <- GenomeInfoDb::seqlevelsStyle(genome)
        ref_organism <- GenomeInfoDb::organism(genome)
        
        if (group == "auto+sex"){
            groups <- c(GenomeInfoDb::extractSeqlevelsByGroup(species = ref_organism,
                                                              style = ref_style,
                                                              group = "auto"),
                        GenomeInfoDb::extractSeqlevelsByGroup(species = ref_organism,
                                                              style = ref_style,
                                                              group = "sex"))
            
            # In some cases, the seqlevelsStyle returns multiple styles.
            # In this case, we need to do a little more work to extract
            # a vector of seqlevels from it.
            groups_names <- names(groups)
            if (! is.null(groups_names))
            {
                # The seqlevels in the groups are now duplicated.
                # The following code deduplicates the list items, so that
                # creating a data frame will work as expected.
                unique_names <- unique(groups_names)
                groups <- plyr::llply(unique_names, function(x) groups[groups_names == x])
                groups <- plyr::llply(groups, unlist, recursive = FALSE)
                
                # In case there are multiple styles applied, we only use the first.
                groups <- unique(as.vector(groups[[1]]))
            }
        }
        else
        {
            groups <- GenomeInfoDb::extractSeqlevelsByGroup ( species = ref_organism,
                                                              style = ref_style,
                                                              group = group )
            groups <- unique(as.vector(t(groups)))
        }
        
        # The provided VCF files may not contain all chromosomes that are
        # available in the reference genome.  Therefore, we only take the
        # chromosomes that are actually available in the VCF file,
        # belonging to the filter group.
        groups <- BiocGenerics::intersect(groups, seqlevels(gr))
        
        # We use 'pruning.mode = "tidy"' to minimize the deleterious effect
        # on variants, yet, remove all variants that aren't in the filter
        # group.  By default, keepSeqlevels would produce an error.
        gr <- GenomeInfoDb::keepSeqlevels(gr, groups, pruning.mode = "tidy")
    }
    return(gr)
}