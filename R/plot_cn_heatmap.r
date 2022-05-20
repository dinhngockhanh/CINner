# =================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
plot_cn_heatmap <- function(model = "",
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

            sample_genotype_profiles <- package_sample$sample_genotype_profiles
            phylogeny_clustering_truth <- package_sample_phylogeny$phylogeny_clustering_truth

            sample_clustering <- phylogeny_clustering_truth$clustering
            sample_tree <- phylogeny_clustering_truth$tree
        }
        #-----------------------------------------Plot total CN profiles
        if (plotcol == "total-copy") {
            filename <- paste(model, "_sim", i, "_CN_total", ".jpeg", sep = "")
            jpeg(file = filename, width = width, height = height)
            p <- plotHeatmap(sample_genotype_profiles,
                plotcol = "state",
                clusters = sample_clustering,
                tree = sample_tree,
                reorderclusters = TRUE,
                plottree = TRUE,
                plotfrequency = TRUE
            )
            print(p)
            dev.off()
        }
        #-----------------------------------------Plot minor CN profiles
        if (plotcol == "minor-copy") {
            filename <- paste(model, "_sim", i, "_CN_minor", ".jpeg", sep = "")
            jpeg(file = filename, width = width, height = height)
            p <- plotHeatmap(sample_genotype_profiles,
                plotcol = "Min",
                clusters = sample_clustering,
                tree = sample_tree,
                reorderclusters = TRUE,
                plottree = TRUE,
                plotfrequency = TRUE
            )
            print(p)
            dev.off()
        }
    }
}
