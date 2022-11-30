# =================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
#' @export
plot_cn_heatmap <- function(model = "",
                            n_simulations = 0,
                            plotcol = "",
                            CN_data = "TRUTH",
                            phylo = "TRUTH",
                            filename_suffix = NULL,
                            width = 1000,
                            height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #-----------------------------------------Decide filename suffix
        if (is.null(filename_suffix)) {
            if (CN_data == "TRUTH") {
                filename_suffix <- "_cnGROUNDTRUTH"
            } else if (CN_data == "NEUTRAL-VARIATIONS") {
                filename_suffix <- "_cnNEUVAR"
            } else if (CN_data == "HMM") {
                filename_suffix <- "_cnGROUNDTRUTH+HMM"
            } else if (CN_data == "NEUTRAL-VARIATIONS-HMM") {
                filename_suffix <- "_cnNEUVAR+HMM"
            }
            if (phylo == "TRUTH") {
                filename_suffix <- paste(filename_suffix, "_phyloGROUNDTRUTH", sep = "")
            } else if (phylo == "UMAP") {
                filename_suffix <- paste(filename_suffix, "_phyloUMAP", sep = "")
            }
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
            if (length(vec_delete) > 0) {
                sample_genotype_profiles <- sample_genotype_profiles[-vec_delete, ]
            }
        } else if (CN_data == "NEUTRAL-VARIATIONS") {
            #---Extract CN profiles from GROUND TRUTH + NEUTRAL VARIATIONS
            sample_genotype_profiles <- simulation$neutral_variations$sample$cn_profiles_long
            #   Remove internal nodes from true CN profiles
            vec_delete <- c(
                grep("Internal-node-", sample_genotype_profiles$cell_id, value = FALSE),
                grep("Initial-clone-", sample_genotype_profiles$cell_id, value = FALSE)
            )
            if (length(vec_delete) > 0) {
                sample_genotype_profiles <- sample_genotype_profiles[-vec_delete, ]
            }
        } else if (CN_data == "HMM") {
            #---Extract CN profiles from GROUND TRUTH + HMMcopy
            sample_genotype_profiles <- simulation$sample$cn_profiles_long_hmm
        } else if (CN_data == "NEUTRAL-VARIATIONS-HMM") {
            #---Extract CN profiles from GROUND TRUTH + NEUTRAL VARIATIONS + HMMcopy
            sample_genotype_profiles <- simulation$neutral_variations$sample$cn_profiles_long_hmm
        }
        #-------------------------------Extract clustering and phylogeny
        if (phylo == "TRUTH") {
            #---Extract clustering and phylogeny from GROUND TRUTH
            sample_clustering <- simulation$sample_phylogeny$phylogeny_clustering_truth$clustering
            sample_tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
        } else if (phylo == "UMAP") {
            #---Find clustering and phylogeny from UMAP
            pctcells <- 0.05
            umapmetric <- "euclidean"
            seed <- NULL
            ncells <- length(unique(sample_genotype_profiles$cell_id))
            CNbins <- sample_genotype_profiles

            clustering_results <- umap_clustering(CNbins,
                minPts = max(round(pctcells * ncells), 2),
                field = "copy",
                umapmetric = umapmetric,
                seed = seed
            )
            sample_clustering <- clustering_results$clustering %>% dplyr::select(cell_id, clone_id)
            sample_tree <- clustering_results$tree
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
