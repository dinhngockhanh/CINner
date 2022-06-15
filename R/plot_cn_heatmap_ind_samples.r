plot_cn_heatmap_ind_samples <- function(model = "",
                                        n_simulations = 0,
                                        plotcol = "",
                                        phylo = TRUE,
                                        width = 1000,
                                        height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #----------Extract CN profiles and true clustering and phylogeny
        if (phylo == TRUE) {
            package_sample <- simulation$sample
            package_sample_phylogeny <- simulation$sample_phylogeny

            sample_genotype_profiles <- package_sample$cn_profiles_long
            phylogeny_clustering_truth <- package_sample_phylogeny$phylogeny_clustering_truth

            sample_clustering <- phylogeny_clustering_truth$clustering
            sample_tree <- phylogeny_clustering_truth$tree
        }
        #-------------------------------------Add column for sample ID's
        sample_genotype_profiles$sample_id <- sub("-.*", "", sample_genotype_profiles$cell_id)
        list_sample_id <- unique(sample_genotype_profiles$sample_id)
        list_cell_id <- unique(sample_genotype_profiles$cell_id)
        #------------------------Split CN profiles into separate samples
        list_sample_genotype_profiles <- list()
        for (j in 1:length(list_sample_id)) {
            sample_id <- list_sample_id[j]
            list_sample_genotype_profiles[[j]] <- sample_genotype_profiles[sample_genotype_profiles$sample_id == sample_id, ]
            list_sample_genotype_profiles[[j]] <- subset(list_sample_genotype_profiles[[j]], select = -c(sample_id))
        }
        #---------------Split phylogeny clustering into separate samples
        sample_clustering$sample_id <- sub("-.*", "", sample_clustering$cell_id)
        list_sample_clustering <- list()
        for (j in 1:length(list_sample_id)) {
            sample_id <- list_sample_id[j]
            list_sample_clustering[[j]] <- sample_clustering[sample_clustering$sample_id == sample_id, ]
            list_sample_clustering[[j]] <- subset(list_sample_clustering[[j]], select = -c(sample_id))
        }
        #---------------------Split phylogeny tree into separate samples
        list_sample_tree <- list()
        for (j in 1:length(list_sample_id)) {
            #   Find cells in this sample
            sample_id <- list_sample_id[j]
            sample_cell_id <- c()
            for (k in 1:length(list_cell_id)) {
                if (sub("-.*", "", list_cell_id[k]) == sample_id) {
                    sample_cell_id <- c(sample_cell_id, list_cell_id[k])
                }
            }
            #   Get subtree for these cells
            list_sample_tree[[j]] <- ape::keep.tip(sample_tree, tip = sample_cell_id)
        }
        #-------------------------Plot total CN profiles for each sample
        if (plotcol == "total-copy") {
            for (j in 1:length(list_sample_id)) {
                filename <- paste(model, "_sim", i, "_CN_total_", list_sample_id[j], ".jpeg", sep = "")
                jpeg(file = filename, width = width, height = height)
                p <- plotHeatmap(list_sample_genotype_profiles[[j]],
                    plotcol = "state",
                    clusters = list_sample_clustering[[j]],
                    tree = list_sample_tree[[j]],
                    reorderclusters = TRUE,
                    plottree = TRUE,
                    plotfrequency = TRUE
                )
                print(p)
                dev.off()
            }
        }
        #-------------------------Plot minor CN profiles for each sample
        if (plotcol == "minor-copy") {
            for (j in 1:length(list_sample_id)) {
                filename <- paste(model, "_sim", i, "_CN_minor_", list_sample_id[j], ".jpeg", sep = "")
                jpeg(file = filename, width = width, height = height)
                p <- plotHeatmap(list_sample_genotype_profiles[[j]],
                    plotcol = "Min",
                    clusters = list_sample_clustering[[j]],
                    tree = list_sample_tree[[j]],
                    reorderclusters = TRUE,
                    plottree = TRUE,
                    plotfrequency = TRUE
                )
                print(p)
                dev.off()
            }
        }
    }
}
