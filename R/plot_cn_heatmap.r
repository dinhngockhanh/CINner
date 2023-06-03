# =================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
#' @export
plot_cn_heatmap <- function(model = "",
                            n_simulations = 0,
                            folder_workplace = "",
                            folder_plots = "",
                            plotcol = "",
                            plottree = TRUE,
                            plotfrequency = TRUE,
                            show_legend = TRUE,
                            show_library_label = TRUE,
                            show_clone_label = TRUE,
                            CN_data = "TRUTH",
                            phylo = "TRUTH",
                            filename_suffix = NULL,
                            width = 1000,
                            height = 500,
                            compute_parallel = TRUE,
                            n_cores = NULL,
                            R_libPaths = NULL) {
    library(signals)
    if (compute_parallel == FALSE) {
        #-----------------------------Plot CN heatmap in sequential mode
        for (iteration in 1:n_simulations) {
            plot_cn_heatmap_one_simulation(
                model,
                iteration,
                folder_workplace,
                folder_plots,
                plotcol,
                plottree,
                plotfrequency,
                show_legend,
                show_library_label,
                show_clone_label,
                CN_data,
                phylo,
                filename_suffix,
                width,
                height
            )
        }
    } else {
        #-------------------------------Plot CN heatmap in parallel mode
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        if (is.null(R_libPaths) == FALSE) {
            R_libPaths <<- R_libPaths
            clusterExport(cl, varlist = c("R_libPaths"))
            clusterEvalQ(cl = cl, .libPaths(R_libPaths))
        }
        clusterEvalQ(cl = cl, library(signals))
        #   Prepare input parameters for plotting
        plot_cn_heatmap_one_simulation <<- plot_cn_heatmap_one_simulation
        model <<- model
        folder_workplace <<- folder_workplace
        folder_plots <<- folder_plots
        plotcol <<- plotcol
        plottree <<- plottree
        plotfrequency <<- plotfrequency
        show_legend <<- show_legend
        show_library_label <<- show_library_label
        show_clone_label <<- show_clone_label
        CN_data <<- CN_data
        phylo <<- phylo
        filename_suffix <<- filename_suffix
        width <<- width
        height <<- height
        clusterExport(cl, varlist = c(
            "plot_cn_heatmap_one_simulation",
            "model",
            "folder_workplace",
            "folder_plots",
            "plotcol",
            "plottree",
            "plotfrequency",
            "show_legend",
            "show_library_label",
            "show_clone_label",
            "CN_data",
            "phylo",
            "filename_suffix",
            "width",
            "height"
        ))
        #   Plot in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
            plot_cn_heatmap_one_simulation(
                model,
                iteration,
                folder_workplace,
                folder_plots,
                plotcol,
                plottree,
                plotfrequency,
                show_legend,
                show_library_label,
                show_clone_label,
                CN_data,
                phylo,
                filename_suffix,
                width,
                height
            )
        })
        #   Stop parallel cluster
        stopCluster(cl)
    }
}

plot_cn_heatmap_one_simulation <- function(model,
                                           iteration,
                                           folder_workplace,
                                           folder_plots,
                                           plotcol,
                                           plottree,
                                           plotfrequency,
                                           show_legend,
                                           show_library_label,
                                           show_clone_label,
                                           CN_data,
                                           phylo,
                                           filename_suffix,
                                           width,
                                           height) {
    #------------------------------------------Input simulation file
    filename <- paste(folder_workplace, model, "_simulation_", iteration, ".rda", sep = "")
    load(filename)
    #-----------------------------------------Decide filename suffix
    if (is.null(filename_suffix)) {
        if (CN_data == "TRUTH") {
            filename_suffix <- "_cnTRUTH"
        } else if (CN_data == "NEUTRAL-VARIATIONS") {
            filename_suffix <- "_cnNEUVAR"
        } else if (CN_data == "HMM") {
            filename_suffix <- "_cnTRUTH+HMM"
        } else if (CN_data == "NEUTRAL-VARIATIONS-HMM") {
            filename_suffix <- "_cnNEUVAR+HMM"
        }
        if (phylo == "TRUTH") {
            filename_suffix <- paste(filename_suffix, "_phyloTRUTH", sep = "")
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
    sample_genotype_profiles$chr <- as.character(sample_genotype_profiles$chr)
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
        filename <- paste(folder_plots, model, "_sim", iteration, "_CN_total", filename_suffix, ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- plotHeatmap(sample_genotype_profiles,
            plotcol = "state",
            clusters = sample_clustering,
            tree = sample_tree,
            reorderclusters = TRUE,
            plottree = plottree,
            plotfrequency = plotfrequency,
            show_legend = show_legend,
            show_library_label = show_library_label,
            show_clone_label = show_clone_label,
            normalize_tree = FALSE,
        )
        print(p)
        dev.off()
    }
    #-----------------------------------------Plot minor CN profiles
    if (plotcol == "minor-copy") {
        filename <- paste(folder_plots, model, "_sim", iteration, "_CN_minor", filename_suffix, ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- plotHeatmap(sample_genotype_profiles,
            plotcol = "Min",
            clusters = sample_clustering,
            tree = sample_tree,
            reorderclusters = TRUE,
            plottree = plottree,
            plotfrequency = plotfrequency,
            show_legend = show_legend,
            show_library_label = show_library_label,
            show_clone_label = show_clone_label,
            normalize_tree = FALSE
        )
        print(p)
        dev.off()
    }
}
