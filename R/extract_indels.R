#' Extract indel signatures
#'
#' @description Will return a 1-column matrix containing the absolute indel signature
#' contributions (i.e. the number of mutations contributing to each mutational signature). The
#' signatures used are insertions/deletions within repeat regions (ins.rep, del.rep),
#' insertions/deletions with flanking microhomology (ins.mh, del.mh), and insertions/deletions
#' which don't fall under the previous 2 categories (ins.none, del.none). Each category is further
#' stratified by the length of the indel.
#'
#' @param bed A dataframe containing the columns: chrom, pos, ref, alt
#' @param sample.name If a character is provided, the header for the output matrix will be named to this. If none is
#' provided, the basename of the vcf file will be used.
#' @param ref_genome A character naming the BSgenome reference genome. Default is "BSgenome.Hsapiens.UCSC.hg19". If another
#' reference genome is indicated, it will also need to be installed.
#' @param indel.len.cap Specifies the max indel sequence length to consider when counting 'repeat' and 'none' contexts.
#' Counts of longer indels will simply be binned to the counts of contexts at the max indel sequence length.
#' @param n.bases.mh.cap Specifies the max bases in microhomology to consider when counting repeat and microhomology
#' contexts. Counts of longer indels will simply be binned to the counts of contexts at the max indel sequence length.
#' @param get.other.indel.allele Only applies when mode=='indel' For indels, some vcfs only report
#' the sequence of one allele (REF for deletions and ALT for insertions). If TRUE, the unreported
#' allele will be retrieved from the genome: a 5' base relative to the indel sequence. This base
#' will also be added to the indel sequence and the POS will be adjusted accordingly (POS=POS-1).
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to variantsFromVcf()
#'
#' @return A 1-column matrix
#' @export
extract_indels <- function(bed, context.database, sample.name=NULL, ref.genome=DEFAULT_GENOME,
                           indel.len.cap=5, n.bases.mh.cap=5, verbose=F, ...)
{
    df <- get_contexts_indel(bed, ...)
    
    if(verbose){ message('Initializing indel signature output vector...') }
    indel_sig_names <- c(
      paste0('del.rep.len.', 1:indel.len.cap),
      paste0('ins.rep.len.', 1:indel.len.cap),
      paste0('del.mh.bimh.', 1:n.bases.mh.cap),
      paste0('ins.mh.bimh.', 1:n.bases.mh.cap),
      paste0('del.none.len.', 1:indel.len.cap),
      paste0('ins.none.len.', 1:indel.len.cap)
    )
    indel_sigs <- structure(rep(0,length(indel_sig_names)), names=indel_sig_names)
    
    # Don't process empty vcfs (df==NA if empty)
    if(!is.data.frame(df)){ stop("Vcf is empty") }
    
    #--------- Pre-calculations for repeat and microhomology contexts ---------#
    if(verbose){ message('Determining the start/end positions for the left/right flanks of each indel...') }
  
    ## Getting flank start/end positions first, then retrieving the sequences using getSeq with a
    ## vectorized input improves speed significantly.
    ## Cap n.indel.lengths.r to 3 (length of 3' (right-hand side) sequence to retrieve),
    ## since the max n_copies_along_flank condition used below caps at >=2.
    ## This improves speed significantly
    flanks_start_end <- with(df,{
      do.call(rbind, Map(
        f = indel_seq_flanks_start_end,
        chrom, pos, indel_len, indel_type,
        n.indel.lengths.r=3
      ))
    })
    
    if(verbose){ message('Retrieving flanking sequences...') }
    l_flank <- getSeq(
      x = eval(parse(text=ref.genome)),
      names = df$chrom,
      start = flanks_start_end[,'l_start'],
      end = flanks_start_end[,'l_end'],
      as.character = T
    )
    
    r_flank <- getSeq(
      x = eval(parse(text=ref.genome)),
      names = df$chrom,
      start = flanks_start_end[,'r_start'],
      end = flanks_start_end[,'r_end'],
      as.character = T
    )
    
    #--------- Repeat contexts ---------#
    if(verbose){ message("Calculating the number of copies of the indel sequence are present in the 3' flanking sequence...") }
    n_copies_along_flank <- unlist(Map(n_copies_along_flank, df$indel_seq, r_flank, USE.NAMES=F))
  
    #--------- Microhomology contexts ---------#
    if(verbose){ message("Calculating the (max) number of bases that are homologous to the 5'/3' flanking sequence...") }
    n_bases_mh <- unlist(Map(function(indel_seq, l_flank, r_flank){
      mh_l <- n_bases_mh(reverse(indel_seq), reverse(l_flank))
      mh_r <- n_bases_mh(indel_seq, r_flank)
      
      max(mh_l,mh_r)
      
    }, df$indel_seq, l_flank, r_flank, USE.NAMES=F))
    
    df <- cbind(df, "repeats"=n_copies_along_flank, "bimh"=n_bases_mh)
    
    # For native indel context 
    if (context.database == "native")
    {
      #--------- Assign repeat, microhomology, or no context ---------#
      if(verbose){ message('Determining indel contexts...') }
      context <- unlist(Map(function(n_copies_along_flank, n_bases_mh, indel_len){
        if (n_copies_along_flank >= 2){
          if(indel_len < 50){ context <-'rep' }
          else { context <- 'mh' }
          
        } else if(n_copies_along_flank >= 1 && n_bases_mh >= 2) {
          context <- 'mh'
        } else if(n_copies_along_flank >= 1 && n_bases_mh >= 1 && indel_len > 3 ) {
          context <- 'mh'
          
        } else {
          context <- 'none'
        }
        
        return(context)
        
      }, n_copies_along_flank, n_bases_mh, df$indel_len))
      
      #--------- Gather components for counting final contexts/signatures ---------#
      if(verbose){ message('Counting indel context occurrences...') }
      ## Slightly redundant (could have just assigned components to a dataframe). But easier to debug
      sig_parts <- data.frame(
        indel_type = df$indel_type,
        context,
        indel_len = df$indel_len,
        n_copies_along_flank,
        n_bases_mh
      )
      
      ## Bin values larger than cap into one bin for indel_len and n_bases_mh
      sig_parts <- within(sig_parts,{
        indel_len[indel_len >= indel.len.cap] <- indel.len.cap
        n_bases_mh[n_bases_mh >= n.bases.mh.cap] <- n.bases.mh.cap
      })
      
      contexts <- with(sig_parts,{
        unlist(Map(function(indel_type,context,indel_len,n_copies_along_flank,n_bases_mh){
          if(context == 'mh'){
            paste(indel_type, context, 'bimh', n_bases_mh, sep = '.')
          } else {
            paste(indel_type, context, 'len', indel_len, sep = '.')
          }
        },
        indel_type, context, indel_len, n_copies_along_flank, n_bases_mh))
      })
      
      return(contexts)
    } else if (context.database == "cosmic")
    {
      # For cosmic indel context
      df$final_context = NA
      
      for (i in 1:nrow(df))
      {
        df$final_context[i] = cosmic_indel_context(df[i,])
      }
      
      df = df[df$final_context %in% indel_context,]
      
      return(df$final_context)
    }
    
    ## Count occurrences of each signature
    # sig_occurrences <- table(with(sig_parts,{
    #   unlist(Map(function(indel_type,context,indel_len,n_copies_along_flank,n_bases_mh){
    #     if(context == 'mh'){
    #       paste(indel_type, context, 'bimh', n_bases_mh, sep = '.')
    #     } else {
    #       paste(indel_type, context, 'len', indel_len, sep = '.')
    #     }
    #   },
    #   indel_type, context, indel_len, n_copies_along_flank, n_bases_mh))
    # }))
    
    ## Fill in indel signature matrix that was initiated at the start of the function
    # indel_sigs[names(sig_occurrences)] <- sig_occurrences
  
    
    # if(verbose){ message('Returning indel context counts...') }
    # out <- matrix(indel_sigs, ncol = 1)
    # rownames(out) <- names(indel_sigs)
    #colnames(out) <- if(is.null(sample.name)){ basename(vcf.file) } else { sample.name }
    
    # return(out)
}
