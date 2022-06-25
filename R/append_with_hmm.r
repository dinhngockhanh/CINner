append_with_hmm <- function(model = "",
                            n_simulations = 0,
                            UMAP = TRUE,
                            pseudo_corrected_readcount = FALSE) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        sample_cell_ID <- simulation$sample$sample_cell_ID
        all_sample_genotype <- simulation$sample$all_sample_genotype
        #-------------Find the pseudo-corrected readcounts for each cell
        #-------------------------inferred by HMMcopy in the wide format
        if (pseudo_corrected_readcount == TRUE) {
            for (cell in 1:length(sample_cell_ID)) {
                #   Find the CN profile for this cell inferred by HMMcopy
                cell_filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", i, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
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
            filename_pseudo_corrected_readcount <- paste(model, "_pseudo_corrected_readcounts_wide_", i, ".csv", sep = "")
            write.csv(pseudo_corrected_readcount_wide_hmm, filename_pseudo_corrected_readcount, row.names = FALSE)
        }
        #-----------------------------Find the CN profiles for each cell
        #-------------------------inferred by HMMcopy in the long format
        hmm_profiles_long_list <- vector("list", length = length(all_sample_genotype))
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", i, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Add record for this cell to list
            hmm_profiles_long_list[[cell]] <- cell_genotype_profile_hmm
        }
        cn_profiles_long_hmm <- rbindlist(hmm_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
        class(cn_profiles_long_hmm) <- "data.frame"
        simulation$sample$cn_profiles_long_hmm <- cn_profiles_long_hmm
        #------------------------Find clustering and phylogeny with UMAP
        if (UMAP == TRUE) {
            pctcells <- 0.05
            umapmetric <- "euclidean"
            seed <- NULL
            ncells <- length(unique(cn_profiles_long_hmm$cell_id))
            CNbins <- cn_profiles_long_hmm

            clustering_results <- umap_clustering(CNbins,
                minPts = max(round(pctcells * ncells), 2),
                field = "copy",
                umapmetric = umapmetric,
                seed = seed
            )
            sample_tree <- clustering_results$tree
            sample_clustering <- clustering_results$clustering %>%
                dplyr::select(cell_id, clone_id)

            phylogeny_clustering_umap_on_hmm <- list()
            phylogeny_clustering_umap_on_hmm$clustering <- sample_clustering
            phylogeny_clustering_umap_on_hmm$tree <- sample_tree
            simulation$sample_phylogeny$phylogeny_clustering_umap_on_hmm <- phylogeny_clustering_umap_on_hmm
        }
        #-----------------------------------------Save simulation record
        save(simulation, file = filename)
    }
}
