#' Fit mutational signatures to a mutation matrix
#' 
#' Find a linear non-negative combination of mutation signatures that
#' reconstructs the mutation matrix. First an optimal reconstruction is achieved via `fit_to_signatures`.
#' However, this is prone to overfitting. 
#' To solve this the signature with the lowest contribution is removed and refitting is repeated.
#' This is done in an iterative fashion. 
#' Each time the cosine similarity between the original and reconstructed profile is calculated.
#' Iterations are stopped when the difference between two iterations becomes more than `max_delta`.
#' The second-last set of signatures is then used for a final refit.
#' 
#' 
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param signatures Signature matrix (dimensions: x mutation types
#' X n signatures)
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity between two iterations. 
#' @return A list containing a fit_res object, similar to `fit_to_signatures` and a list of ggplot graphs
#' that for each sample shows in what order the signatures were removed and how this affected the cosine similarity.
#' 
#' @seealso \code{\link{mut_matrix}}, 
#' \code{\link{fit_to_signatures}}, 
#' \code{\link{get_cos_sim_ori_vs_rec}}, 
#' \code{\link{plot_sim_decay}} 
#' @export
#' 
#' @importFrom magrittr %>% 
#'
#' @examples
#' #' ## See the 'mut_matrix()' example for how we obtained the mutation matrix:
#' mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#'                     package="MutationalPatterns"))
#'
#' ## You can download the signatures from the COSMIC website:
#' # https://cancer.sanger.ac.uk/cosmic/signatures
#' 
#' ## We copied the file into our package for your convenience.
#' filename <- system.file("extdata/snv_signatures_probabilities.txt",
#'                         package="MutationalPatterns")
#' signatures <- read.table(filename, sep = "\t", header = TRUE)
#' 
#' ## Remove unnecessary columns
#' signatures = as.matrix(signatures[,-c(1,2)])
#' 
#' ## Fit to signatures strict
#' strict_refit = fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.05)
#' 
#' ## fit_res similar to `fit_to_signatures`
#' fit_res = strict_refit$fit_res
#' 
#' ## list of ggplots that shows how the cosine similarity was reduced during the iterations
#' fig_l = strict_refit$sim_decay_fig
fit_to_signatures_strict = function(mut_matrix, signatures, max_delta = 0.05){
    
    #Remove signatures with zero contribution across samples
    fit_res = fit_to_signatures(mut_matrix, signatures)
    sig_pres = rowSums(fit_res$contribution) != 0
    my_signatures_total = signatures[,sig_pres, drop = F]
    nsigs = ncol(my_signatures_total)
    
    #perform signature selection per sample
    all_results = vector("list", ncol(mut_matrix))
    for (i in seq(1, ncol(mut_matrix))){
        my_signatures = my_signatures_total
        mut_mat_sample = mut_matrix[, i, drop = F]
        
        #Fit again
        fit_res = fit_to_signatures(mut_mat_sample, my_signatures)
        sim = get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
        
        #Keep track of the cosine similarity and which signatures are removed.
        sims = vector("list", nsigs)
        sims[[1]] = sim
        removed_sigs = vector("list", nsigs)
        removed_sigs[[1]] = "None"
        
        #Sequentially remove the signature with the lowest contribution
        for (j in seq(2, nsigs)){
            
            #Remove signature with the weakest relative contribution
            contri_order = fit_res$contribution %>% 
                prop.table(2) %>% 
                rowSums() %>% 
                order()
            weakest_sig_index = contri_order[1]
            weakest_sig = colnames(my_signatures)[weakest_sig_index]
            removed_sigs[[j]] = weakest_sig
            signatures_sel = my_signatures[,-weakest_sig_index, drop = F]
            
            
            #Fit with new signature selection
            fit_res = fit_to_signatures(mut_mat_sample, signatures_sel)
            sim_new = get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
            
            if (is.nan(sim_new) == TRUE){
                sim_new = 0
                warning("New similarity between the original and the reconstructed 
                        spectra after the removal of a signature was NaN. 
                        It has been converted into a 0. 
                        This happened with the following fit_res:")
                print(fit_res)
            }
            sims[[j]] = sim_new
            
            #Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
            delta = sim - sim_new
            if(delta <= max_delta){
                my_signatures = signatures_sel
                sim = sim_new
            }
            else{
                break
            }
        }
        
        #Plot how the cosine similarities decayed
        sim_decay_fig = plot_sim_decay(sims, removed_sigs, max_delta)
        
        #Perform final fit on selected signatures
        fit_res = fit_to_signatures(mut_mat_sample, my_signatures) 
        
        #Add data of sample to list.
        results = list("sim_decay_fig" = sim_decay_fig, "fit_res" = fit_res)
        all_results[[i]] = results
    }
    
    #Get decay figs and fit_res in separate lists
    decay_figs = purrr::map(all_results, "sim_decay_fig")
    fit_res = purrr::map(all_results, "fit_res")
    
    #Combine the contribution of all samples and put it in a matrix
    contribution = purrr::map(fit_res, "contribution") %>% 
        purrr::map(function(x) tibble::rownames_to_column(as.data.frame(x))) %>% 
        purrr::reduce(dplyr::full_join, by = "rowname")
    rownames(contribution) = contribution$rowname
    contribution = contribution %>% 
        dplyr::select(-rowname) %>% 
        as.matrix()
    
    #Combine the reconstructed of all samples
    reconstructed = purrr::map(fit_res, "reconstructed") %>% 
        do.call(cbind, .)
    
    #Combine all and return
    fit_res = list("contribution" = contribution, "reconstructed" = reconstructed)
    results = list("sim_decay_fig" = decay_figs, "fit_res" = fit_res)
    return(results)
}



#' Get the cosine similarity between a reconstructed mutation matrix and the original
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param fit_res Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @return Cosine similarity
#'
get_cos_sim_ori_vs_rec = function(mut_mat, fit_res){
    cos_sim_all = cos_sim_matrix(mut_mat, fit_res$reconstructed)
    cos_sim = diag(cos_sim_all)
    mean_cos_sim = mean(cos_sim)
}


#' Plot decay in cosine similarity as signatures are removed.
#'
#'This function is called by fit_to_signatures_strict
#'
#' @param sims List of cosine similarities
#' @param removed_sigs List of iteratively removed signatures
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity.
#'
#' @import ggplot2
#' @importFrom magrittr %>% 
#'
#' @return ggplot object
#'
plot_sim_decay = function(sims, removed_sigs, max_delta){
    
    #
    sims = sims[!S4Vectors::isEmpty(sims)] %>%
        unlist()
    removed_sigs = removed_sigs[!S4Vectors::isEmpty(removed_sigs)] %>% 
        unlist()
    tb = tibble::tibble("Cosine_similarity" = sims, 
                        "Removed_signatures" = factor(removed_sigs, levels = removed_sigs))
    
    #Determine if the final removed signature exceeded the cutoff.
    sims_l = length(sims)
    col = rep("low_delta", sims_l)
    final_delta = sims[sims_l-1] - sims[sims_l]
    if (final_delta > max_delta){
        col[sims_l] = "high_delta"
    }
    
    fig = ggplot(data = tb, aes(x = Removed_signatures, y = Cosine_similarity, fill = col)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(limits = c("low_delta", "high_delta"), values = c("grey", "red"), guide = F) +
        labs(x = "Removed signatures", y = paste0("Cosine similarity (max delta: ", max_delta, ")")) +
        theme(axis.text.x = element_text(angle = 90, size = 12), text = element_text(size=20))
    return(fig)
}
