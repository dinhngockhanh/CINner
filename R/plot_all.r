#' @export
plot_all <- function(model = "",
                     n_simulations = 0,
                     unit_time = "year",
                     folder_workplace = NULL,
                     folder_plots = NULL,
                     compute_parallel = TRUE,
                     n_cores = NULL,
                     R_libPaths = NULL) {
    if (is.null(folder_workplace)) {
        folder_workplace <- ""
    } else {
        dir.create(folder_workplace)
        folder_workplace <- paste(folder_workplace, "/", sep = "")
    }
    if (is.null(folder_plots)) {
        folder_plots <- ""
    } else {
        dir.create(folder_plots)
        folder_plots <- paste(folder_plots, "/", sep = "")
    }
    # #--------------------------------------Plot average ploidy over time
    # cat("Plotting average ploidy over time...\n")
    # start_time <- Sys.time()
    # plot_average_ploidy(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots = folder_plots,
    #     unit_time = unit_time,
    #     width = 2000,
    #     height = 1000,
    #     compute_parallel = compute_parallel,
    #     n_cores = n_cores,
    #     R_libPaths = R_libPaths
    # )
    # end_time <- Sys.time()
    # print(end_time - start_time)
    # cat("\n")



    #---------------------------------Plot clonal evolution as fish plot
    cat("Plotting clonal evolution as fish plot...\n")
    start_time <- Sys.time()
    plot_clonal_fishplot(
        model = model,
        n_simulations = n_simulations,
        folder_workplace = folder_workplace,
        folder_plots = folder_plots,
        unit_time = unit_time,
        width = 2000,
        height = 1000,
        compute_parallel = compute_parallel,
        n_cores = n_cores,
        R_libPaths = R_libPaths
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    cat("\n")
    #----------------------------Plot clonal evolution as phylogeny tree
    cat("Plotting clonal evolution as phylogeny tree...\n")
    start_time <- Sys.time()
    plot_clonal_phylo(
        model = model,
        n_simulations = n_simulations,
        folder_workplace = folder_workplace,
        folder_plots = folder_plots,
        width = 2000,
        height = 2000,
        compute_parallel = compute_parallel,
        n_cores = n_cores,
        R_libPaths = R_libPaths
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    cat("\n")
    #-----Plot total CN profile - CN profiles = TRUTH, phylogeny = TRUTH
    cat("Plotting total CN profile - CN profiles = TRUTH, phylogeny = TRUTH...\n")
    start_time <- Sys.time()
    plot_cn_heatmap(
        model = model,
        n_simulations = n_simulations,
        folder_workplace = folder_workplace,
        folder_plots = folder_plots,
        plotcol = "total-copy",
        CN_data = "TRUTH",
        phylo = "TRUTH",
        width = 1000,
        height = 1000,
        compute_parallel = compute_parallel,
        n_cores = n_cores,
        R_libPaths = R_libPaths
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    cat("\n")
    # #----Plot total CN profile - CN profiles = NEUVAR, phylogeny = TRUTH
    # cat("Plotting total CN profile - CN profiles = NEUVAR, phylogeny = TRUTH...\n")
    # start_time <- Sys.time()
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots = folder_plots,
    #     plotcol = "total-copy",
    #     CN_data = "NEUTRAL-VARIATIONS",
    #     phylo = "TRUTH",
    #     width = 1000,
    #     height = 1000,
    #     compute_parallel = compute_parallel,
    #     n_cores = n_cores,
    #     R_libPaths = R_libPaths
    # )
    # end_time <- Sys.time()
    # print(end_time - start_time)
    # cat("\n")
    # #-----Plot total CN profile - CN profiles = NEUVAR, phylogeny = UMAP
    # cat("Plotting total CN profile - CN profiles = NEUVAR, phylogeny = UMAP...\n")
    # start_time <- Sys.time()
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots = folder_plots,
    #     plotcol = "total-copy",
    #     CN_data = "NEUTRAL-VARIATIONS",
    #     phylo = "UMAP",
    #     width = 1000,
    #     height = 1000,
    #     compute_parallel = compute_parallel,
    #     n_cores = n_cores,
    #     R_libPaths = R_libPaths
    # )
    # end_time <- Sys.time()
    # print(end_time - start_time)
    # cat("\n")







    # #-------------------------------Plot minor CN profile - GROUND TRUTH
    # cat("Plotting minor CN profile - GROUND TRUTH...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots = folder_plots,
    #     plotcol = "minor-copy",
    #     CN_data = "TRUTH",
    #     phylo = "TRUTH",
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths = R_libPaths
    # )
    # #-------------Plot total CN profile - CN from HMM & phylo from TRUTH
    # cat("Plotting total CN profile - CN=HMMcopy, phylogeny=TRUTH...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots=folder_plots,
    #     plotcol = "total-copy",
    #     CN_data = "HMM",
    #     phylo = "TRUTH",
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths=R_libPaths
    # )
    # #--------------Plot total CN profile - CN from HMM & phylo from UMAP
    # cat("Plotting total CN profile - CN=HMMcopy, phylogeny=UMAP...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots=folder_plots,
    #     plotcol = "total-copy",
    #     CN_data = "HMM",
    #     phylo = "UMAP",
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths=R_libPaths
    # )







    # #----------------------Plot CN profile for each clone - GROUND TRUTH
    # plot_cn_per_clone(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots=folder_plots,
    #     CN_data = "TRUTH",
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths=R_libPaths
    # )
    # #-------------------Plot total CN profile for for individual samples
    # cat("Plotting total CN profile for each sample...\n")
    # plot_cn_heatmap_ind_samples(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots=folder_plots,
    #     plotcol = "total-copy",
    #     phylo = TRUE,
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths=R_libPaths
    # )
    # #-------------------Plot minor CN profile for for individual samples
    # cat("Plotting minor CN profile for each sample...\n")
    # plot_cn_heatmap_ind_samples(
    #     model = model,
    #     n_simulations = n_simulations,
    #     folder_workplace = folder_workplace,
    #     folder_plots=folder_plots,
    #     plotcol = "minor-copy",
    #     phylo = TRUE,
    #     width = 1000,
    #     height = 1000,
    #     R_libPaths=R_libPaths
    # )
    #----------------------------------------------Plot total population
    cat("Plotting total population vs input dynamics...\n")
    start_time <- Sys.time()
    plot_tot_pop_vs_input(
        model = model,
        n_simulations = n_simulations,
        folder_workplace = folder_workplace,
        folder_plots = folder_plots,
        unit_time = unit_time,
        width = 1000,
        height = 500,
        compute_parallel = compute_parallel,
        n_cores = n_cores,
        R_libPaths = R_libPaths
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    cat("\n")
}
