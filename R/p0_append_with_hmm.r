#' @export
p0_append_with_hmm <- function(simulation,
                               model = "",
                               iteration = 0,
                               pseudo_corrected_readcount = FALSE,
                               model_readcount_base = "all") {
    sample_cell_ID <- simulation$sample$sample_cell_ID
    #----------------Find the pseudo-corrected readcounts in wide format
    #-----------------------inferred by HMMcopy from true-CN-based reads
    if ((pseudo_corrected_readcount == TRUE) & ((model_readcount_base == "truth") | (model_readcount_base == "all"))) {
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", iteration, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Initiate the readcount profiles with genome location indices
            if (cell == 1) {
                pseudo_corrected_readcount_wide_hmm <- cell_genotype_profile_hmm[, 1:3]
            }
            #   Add column for cell ID
            cell_id <- sample_cell_ID[cell]
            pseudo_corrected_readcount_wide_hmm[cell_id] <- round(mean(cell_genotype_profile_hmm$reads) * cell_genotype_profile_hmm$cor_gc, 0)
        }
        simulation$sample$pseudo_corrected_readcount_wide_hmm <- pseudo_corrected_readcount_wide_hmm
        filename_pseudo_corrected_readcount <- paste(model, "_pseudo_corrected_readcounts_wide_", iteration, ".csv", sep = "")
        write.csv(pseudo_corrected_readcount_wide_hmm, filename_pseudo_corrected_readcount, row.names = FALSE)
    }
    #----------------Find the pseudo-corrected readcounts in wide format
    #---------------------inferred by HMMcopy from neuvar-CN-based reads
    if ((pseudo_corrected_readcount == TRUE) & ((model_readcount_base == "neuvar") | (model_readcount_base == "all"))) {
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_neuvar_cn_profiles_long_", iteration, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Initiate the readcount profiles with genome location indices
            if (cell == 1) {
                pseudo_corrected_readcount_wide_hmm <- cell_genotype_profile_hmm[, 1:3]
            }
            #   Add column for cell ID
            cell_id <- sample_cell_ID[cell]
            pseudo_corrected_readcount_wide_hmm[cell_id] <- round(mean(cell_genotype_profile_hmm$reads) * cell_genotype_profile_hmm$cor_gc, 0)
        }
        simulation$neutral_variations$sample$pseudo_corrected_readcount_wide_hmm <- pseudo_corrected_readcount_wide_hmm
        filename_pseudo_corrected_readcount <- paste(model, "_pseudo_corrected_readcounts_wide_", iteration, ".csv", sep = "")
        write.csv(pseudo_corrected_readcount_wide_hmm, filename_pseudo_corrected_readcount, row.names = FALSE)
    }
    #------------------Find the CN profiles for each cell in long format
    #-----------------------inferred by HMMcopy from true-CN-based reads
    if ((model_readcount_base == "truth") | (model_readcount_base == "all")) {
        hmm_profiles_long_list <- vector("list", length = length(sample_cell_ID))
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", iteration, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Add record for this cell to list
            hmm_profiles_long_list[[cell]] <- cell_genotype_profile_hmm
        }
        cn_profiles_long_hmm <- rbindlist(hmm_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
        class(cn_profiles_long_hmm) <- "data.frame"
        simulation$sample$cn_profiles_long_hmm <- cn_profiles_long_hmm
    }
    #------------------Find the CN profiles for each cell in long format
    #---------------------inferred by HMMcopy from neuvar-CN-based reads
    if ((model_readcount_base == "neuvar") | (model_readcount_base == "all")) {
        hmm_profiles_long_list <- vector("list", length = length(sample_cell_ID))
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_neuvar_cn_profiles_long_", iteration, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Add record for this cell to list
            hmm_profiles_long_list[[cell]] <- cell_genotype_profile_hmm
        }
        cn_profiles_long_hmm <- rbindlist(hmm_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
        class(cn_profiles_long_hmm) <- "data.frame"
        simulation$neutral_variations$sample$cn_profiles_long_hmm <- cn_profiles_long_hmm
    }
    return(simulation)
}
