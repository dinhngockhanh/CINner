plot_all <- function(model = "",
                     n_simulations = 0) {
    #   Plot clonal evolution as phylogeny tree
    cat("Plotting clonal evolution as phylogeny tree...\n")
    plot_clonal_phylo(
        model = model,
        n_simulations = n_simulations,
        width = 1000,
        height = 1000
    )
    print("HERE")
    #   Plot cell evolution as phylogeny tree
    cat("Plotting cell evolution as phylogeny tree...\n")
    plot_cell_phylo(
        model = model,
        n_simulations = n_simulations,
        width = 1000,
        height = 1000
    )
    #   Plot total population
    cat("Plotting total population vs input dynamics...\n")
    plot_tot_pop_vs_input(
        model = model,
        n_simulations = n_simulations,
        vec_time = seq(0, 80, by = 1),
        unit_time = "year",
        width = 1000,
        height = 500
    )
    #   Plot total CN profile
    cat("Plotting total CN profile...\n")
    plot_cn_heatmap(
        model = model,
        n_simulations = n_simulations,
        plotcol = "total-copy",
        phylo = TRUE,
        width = 1000,
        height = 1000
    )
    #   Plot minor CN profile
    cat("Plotting minor CN profile...\n")
    plot_cn_heatmap(
        model = model,
        n_simulations = n_simulations,
        plotcol = "minor-copy",
        phylo = TRUE,
        width = 1000,
        height = 1000
    )
    # #   Plot clonal evolution as fish plot
    # cat("Plotting clonal evolution as fish plot...\n")
    # plot_clonal_fishplot(
    #     model = model,
    #     n_simulations = n_simulations,
    #     vec_time = seq(0, 80, by = 1),
    #     unit_time = "year",
    #     width = 1000,
    #     height = 500
    # )
}
