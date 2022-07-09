# =================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
#' @export
plot_cn_heatmap <- function(model = "",
                            n_simulations = 0,
                            plotcol = "",
                            CN_data = "TRUTH",
                            phylo = "TRUTH",
                            width = 1000,
                            height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #-----------------------------------------Decide filename suffix
        if (CN_data == "TRUTH" & phylo == "TRUTH") {
            filename_suffix <- "_cnTRUTH_phyloTRUTH"
        }
        if (CN_data == "HMM" & phylo == "TRUTH") {
            filename_suffix <- "_cnHMM_phyloTRUTH"
        }
        if (CN_data == "HMM" & phylo == "UMAP") {
            filename_suffix <- "_cnHMM_phyloUMAP"
        }
        #--------------------------------------------Extract CN profiles
        if (CN_data == "TRUTH") {
            #---Extract CN profiles from GROUND TRUTH
            sample_genotype_profiles <- simulation$sample$cn_profiles_long
            #   Remove internal nodes from true CN profiles
            vec_delete <- c(
                grep("Internal-node-", sample_genotype_profiles$cell_id, value = FALSE),
                grep("Initial-clone-", sample_genotype_profiles$cell_id, value = FALSE)
            )
            sample_genotype_profiles <- sample_genotype_profiles[-vec_delete, ]
        }
        if (CN_data == "HMM") {
            package_sample <- simulation$sample
            #---Extract CN profiles from HMMcopy
            sample_genotype_profiles <- simulation$sample$cn_profiles_long_hmm
        }
        print(sample_genotype_profiles)
        #-------------------------------Extract clustering and phylogeny
        if (phylo == "TRUTH") {
            #---Extract clustering and phylogeny from GROUND TRUTH
            sample_clustering <- simulation$sample_phylogeny$phylogeny_clustering_truth$clustering
            sample_tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
        }
        if (phylo == "UMAP") {
            #---Find clustering and phylogeny from UMAP
            sample_clustering <- simulation$sample_phylogeny$phylogeny_clustering_umap_on_hmm$clustering
            sample_tree <- simulation$sample_phylogeny$phylogeny_clustering_umap_on_hmm$tree
        }
        #-----------------------------------------Plot total CN profiles
        if (plotcol == "total-copy") {
            filename <- paste(model, "_sim", i, "_CN_total", filename_suffix, ".jpeg", sep = "")
            jpeg(file = filename, width = width, height = height)
            p <- plotHeatmap(sample_genotype_profiles,
                plotcol = "state",
                clusters = sample_clustering,
                tree = sample_tree,
                reorderclusters = TRUE,
                plottree = TRUE,
                normalize_tree = FALSE,
                plotfrequency = TRUE
            )
            print(p)
            dev.off()
        }
        #-----------------------------------------Plot minor CN profiles
        if (plotcol == "minor-copy") {
            filename <- paste(model, "_sim", i, "_CN_minor", filename_suffix, ".jpeg", sep = "")
            jpeg(file = filename, width = width, height = height)
            p <- plotHeatmap(sample_genotype_profiles,
                plotcol = "Min",
                clusters = sample_clustering,
                tree = sample_tree,
                reorderclusters = TRUE,
                plottree = TRUE,
                normalize_tree = FALSE,
            )
            print(p)
            dev.off()
        }
    }
}
