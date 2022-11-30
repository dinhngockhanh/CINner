#' @export
plot_all <- function(model = "",
                     n_simulations = 0) {
    #----------------------------Plot clonal evolution as phylogeny tree
    cat("Plotting clonal evolution as phylogeny tree...\n")
    plot_clonal_phylo(
        model = model,
        n_simulations = n_simulations,
        width = 2000,
        height = 2000
    )
    #----------------------------------------------Plot total population
    cat("Plotting total population vs input dynamics...\n")
    plot_tot_pop_vs_input(
        model = model,
        n_simulations = n_simulations,
        unit_time = "year",
        width = 1000,
        height = 500
    )
    #---------------------------------Plot clonal evolution as fish plot
    cat("Plotting clonal evolution as fish plot...\n")
    plot_clonal_fishplot(
        model = model,
        n_simulations = n_simulations,
        unit_time = "year",
        width = 2000,
        height = 1000
    )



    #-------------------------------Plot total CN profile - ground truth
    cat("Plotting total CN profile - GROUND TRUTH...\n")
    plot_cn_heatmap(
        model = model,
        n_simulations = n_simulations,
        plotcol = "total-copy",
        CN_data = "TRUTH",
        phylo = "TRUTH",
        width = 1000,
        height = 1000
    )
    # #-------------------------------Plot minor CN profile - GROUND TRUTH
    # cat("Plotting minor CN profile - GROUND TRUTH...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     plotcol = "minor-copy",
    #     CN_data = "TRUTH",
    #     phylo = "TRUTH",
    #     width = 1000,
    #     height = 1000
    # )
    # #-------------Plot total CN profile - CN from HMM & phylo from TRUTH
    # cat("Plotting total CN profile - CN=HMMcopy, phylogeny=TRUTH...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     plotcol = "total-copy",
    #     CN_data = "HMM",
    #     phylo = "TRUTH",
    #     width = 1000,
    #     height = 1000
    # )
    # #--------------Plot total CN profile - CN from HMM & phylo from UMAP
    # cat("Plotting total CN profile - CN=HMMcopy, phylogeny=UMAP...\n")
    # plot_cn_heatmap(
    #     model = model,
    #     n_simulations = n_simulations,
    #     plotcol = "total-copy",
    #     CN_data = "HMM",
    #     phylo = "UMAP",
    #     width = 1000,
    #     height = 1000
    # )
    # #----------------------Plot CN profile for each clone - GROUND TRUTH
    # plot_cn_per_clone(
    #     model = model,
    #     n_simulations = n_simulations,
    #     CN_data = "TRUTH",
    #     width = 1000,
    #     height = 1000
    # )
    # #-------------------Plot total CN profile for for individual samples
    # cat("Plotting total CN profile for each sample...\n")
    # plot_cn_heatmap_ind_samples(
    #     model = model,
    #     n_simulations = n_simulations,
    #     plotcol = "total-copy",
    #     phylo = TRUE,
    #     width = 1000,
    #     height = 1000
    # )
    # #-------------------Plot minor CN profile for for individual samples
    # cat("Plotting minor CN profile for each sample...\n")
    # plot_cn_heatmap_ind_samples(
    #     model = model,
    #     n_simulations = n_simulations,
    #     plotcol = "minor-copy",
    #     phylo = TRUE,
    #     width = 1000,
    #     height = 1000
    # )
}
