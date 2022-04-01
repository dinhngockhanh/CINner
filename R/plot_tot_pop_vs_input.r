# ==============PLOT TOTAL POPULATION SIZE - SIMULATION VS INPUT DYNAMICS
plot_tot_pop_vs_input <- function(model = "",
                                  n_simulations = 0,
                                  vec_time = c(0),
                                  unit_time = "year",
                                  width = 1000,
                                  height = 500){
    for (i in 1:n_simulations){
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #----------------Find total population dynamics from input model
        if (unit_time == "day") {
            vec_time <- vec_time
        } else {
            if (unit_time == "week") {
                vec_time <- 7 * vec_time
            } else {
                if (unit_time == "month") {
                    vec_time <- 30 * vec_time
                } else {
                    if (unit_time == "year") {
                        vec_time <- 365 * vec_time
                    }
                }
            }
        }
        SIMULATOR_VARIABLES_for_simulation(model)
        vec_cell_count_exp_plot <- func_expected_population(vec_time)
        #-----------------Find total population dynamics from simulation
        #   Find clonal population records from simulation
        vec_time_simulation <- simulation$clonal_evolution$evolution_traj_time
        vec_clonal_populations <- simulation$clonal_evolution$evolution_traj_population
        #   Find total population at desired time points
        vec_cell_count_sim_plot <- rep(0, length(vec_time))
        for (j in 1:length(vec_time)) {
            time_plot <- vec_time[j]
            loc <- which.min(abs(vec_time_simulation - time_plot))
            vec_cell_count_sim_plot[j] <- sum(vec_clonal_populations[[loc]])
        }
        #-------Plot total population dynamics - simulation vs input expectation
        if (unit_time == "day") {
            vec_time <- vec_time
        } else {
            if (unit_time == "week") {
                vec_time <- vec_time / 7
            } else {
                if (unit_time == "month") {
                    vec_time <- vec_time / 30
                } else {
                    if (unit_time == "year") {
                        vec_time <- vec_time / 365
                    }
                }
            }
        }
        df <- data.frame(vec_time, vec_cell_count_sim_plot, vec_cell_count_exp_plot)

        filename <- paste(model, "_tot_pop_", i, ".jpeg", sep = "")

        jpeg(file = filename, width = width, height = height)

        p <- ggplot(df, aes(vec_time)) +
                geom_line(aes(y = vec_cell_count_exp_plot), colour = "red") +
                geom_point(aes(y = vec_cell_count_sim_plot), colour = "blue") +
                xlim(vec_time[1], vec_time[length(vec_time)]) +
                xlab("Age") +
                ylab("Cell count") +
                ggtitle("Expected (red line) vs Simulated (blue dots) dynamics of total cell count")
        print(p)
        dev.off()
    }
}
