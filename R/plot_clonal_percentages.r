# ===================================PLOT CLONAL EVOLUTION AS PERCENTAGES
# ==================================================AT SAMPLE TIME POINTS
plot_clonal_percentages <- function(model = "",
                                    n_simulations = 0,
                                    width = 1000,
                                    height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        # Input the
    }
}
