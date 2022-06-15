append_with_hmm <- function(model = "",
                            n_simulations = 0,
                            UMAP = TRUE) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        sample_cell_ID <- simulation$sample$sample_cell_ID
        all_sample_genotype <- simulation$sample$all_sample_genotype
        #---------------Find the CN profiles for each cell in the sample
        #--------------------------------------------inferred by HMMcopy
        hmm_profiles_long_list <- vector("list", length = length(all_sample_genotype))
        for (cell in 1:length(sample_cell_ID)) {
            #   Find the CN profile for this cell inferred by HMMcopy
            cell_filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", i, "_", sample_cell_ID[cell], "_hmm_reads.csv", sep = "")
            cell_genotype_profile_hmm <- read.csv(cell_filename)
            #   Add record for this cell to list
            hmm_profiles_long_list[[cell]] <- cell_genotype_profile_hmm
        }
        #---------------------------Bind all cells' CN profiles together
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
