#' @export
plot_phylo_stats <- function(copynumber_sims,
                             stats_DATA = list()) {
    #---------------------Get cell phylogeny statistics from simulations
    treeList <- list()
    for (iteration in 1:n_simulations) {
        treeList[[iteration]] <- copynumber_sims[[iteration]]$cell_phylogeny
    }
    stats_SIMS <- suppressWarnings(phyloTop(treeList))
    #--------------------------------------Create dataframe for plotting
    df_plot <- data.frame(value = c(), statistics = c(), group = c())
    vec_variables <- colnames(stats_SIMS)
    for (i in 1:length(vec_variables)) {
        if (vec_variables[i] %in% colnames(stats_DATA)) {
            variable <- vec_variables[i]

            df_plot_next <- data.frame(
                value = c(
                    stats_SIMS[[variable]],
                    stats_DATA[[variable]]
                ),
                statistics = rep(variable, length = (nrow(stats_SIMS) + nrow(stats_DATA))),
                group = c(
                    rep("simulations", length = nrow(stats_SIMS)),
                    rep("data", length = nrow(stats_DATA))
                )
            )
            df_plot <- rbind(df_plot, df_plot_next)
        }
    }
    #-------------------------------------Create plot for all statistics
    for (i in 1:length(vec_variables)) {
        # for (i in 1:1) {
        variable <- vec_variables[i]
        df_plot_mini <- df_plot[which(df_plot$statistics == variable), ]
        p <- ggplot(df_plot_mini, aes(x = value, fill = group)) +
            geom_density(alpha = .5) +
            xlab(variable) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20))
        print(p)
    }
}
