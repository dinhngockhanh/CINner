# ==============PLOT TOTAL POPULATION SIZE - SIMULATION VS INPUT DYNAMICS
#' @export
plot_tot_pop_vs_input <- function(model = "",
                                  n_simulations = 0,
                                  vec_time = NULL,
                                  unit_time = "year",
                                  width = 1000,
                                  height = 500) {
    library(ggplot2)
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #-----------------------------Transform time points if necessary
        if (is.null(vec_time)) {
            evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
            T_0 <- evolution_traj_time[1]
            T_end <- evolution_traj_time[length(evolution_traj_time)]
            vec_time_expectation <- seq(T_0, T_end, by = ((T_end - T_0) / 100))
            if (unit_time == "day") {
                vec_time <- vec_time_expectation
            } else if (unit_time == "week") {
                vec_time <- vec_time_expectation / 7
            } else if (unit_time == "month") {
                vec_time <- vec_time_expectation / 30
            } else if (unit_time == "year") {
                vec_time <- vec_time_expectation / 365
            }
        } else {
            if (unit_time == "day") {
                vec_time_expectation <- vec_time
            } else if (unit_time == "week") {
                vec_time_expectation <- 7 * vec_time
            } else if (unit_time == "month") {
                vec_time_expectation <- 30 * vec_time
            } else if (unit_time == "year") {
                vec_time_expectation <- 365 * vec_time
            }
        }
        #----------------Find total population dynamics from input model
        # vec_time_expectation <- vec_time
        SIMULATOR_VARIABLES_for_simulation(model)
        vec_cell_count_exp_plot <- func_expected_population(vec_time_expectation)
        #-----------------Find total population dynamics from simulation
        #   Find clonal population records from simulation
        vec_time_simulation <- simulation$clonal_evolution$evolution_traj_time
        vec_clonal_populations <- simulation$clonal_evolution$evolution_traj_population
        #   Find total population at desired time points
        vec_cell_count_sim_plot <- rep(0, length(vec_time))
        for (j in 1:length(vec_time)) {
            time_plot <- vec_time_expectation[j]
            loc <- which.min(abs(vec_time_simulation - time_plot))
            vec_cell_count_sim_plot[j] <- sum(vec_clonal_populations[[loc]])
        }
        #-------Plot total population dynamics - simulation vs input expectation
        df <- data.frame(vec_time, vec_cell_count_sim_plot, vec_cell_count_exp_plot)
        filename <- paste(model, "_sim", i, "_tot_pop", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- ggplot(df, aes(vec_time)) +
            geom_line(aes(y = vec_cell_count_exp_plot), colour = "red", size = 2) +
            geom_point(aes(y = vec_cell_count_sim_plot), colour = "blue", size = 4) +
            xlim(vec_time[1], vec_time[length(vec_time)]) +
            xlab("Age") +
            ylab("Cell count") +
            # ggtitle("Expected (red line) vs Simulated (blue dots) dynamics of total cell count") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20))
        print(p)
        dev.off()
    }
}
