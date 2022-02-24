#==============PLOT TOTAL POPULATION SIZE - SIMULATION VS INPUT DYNAMICS
PLOT_total_population_vs_input <- function(package_simulation,MODEL,vec_time_plot,unit){
    if (unit=='year'){
        vec_time_plot               <- 365*vec_time_plot
    }
#------------------------Find total population dynamics from input model
    SIMULATOR_VARIABLES_for_simulation(MODEL)
    vec_cell_count_exp_plot         <- func_expected_population(vec_time_plot)
#-------------------------Find total population dynamics from simulation
#   Find clonal population records from simulation
    package_clonal_evolution        <- package_simulation[[1]]
    vec_time                        <- package_clonal_evolution[[13]]
    vec_clonal_populations          <- package_clonal_evolution[[16]]
#   Find total population at desired time points
    vec_cell_count_sim_plot         <- rep(0,length(vec_time_plot))
    for (i in 1:length(vec_time_plot)){
        time_plot                   <- vec_time_plot[i]
        loc                         <- which.min(abs(vec_time-time_plot))
        vec_cell_count_sim_plot[i]  <- sum(vec_clonal_populations[[loc]])
    }
#-------Plot total population dynamics - simulation vs input expectation
    if (unit=='year'){
        vec_time_plot               <- vec_time_plot/365
    }
    df                              <- data.frame(vec_time_plot,vec_cell_count_sim_plot,vec_cell_count_exp_plot)
    ggplot(df,aes(vec_time_plot)) +
          geom_line(aes(y=vec_cell_count_exp_plot),colour="red") +
          geom_point(aes(y=vec_cell_count_sim_plot),colour="blue") +
          xlim(vec_time_plot[1],vec_time_plot[length(vec_time_plot)]) +
          xlab('Age') +
          ylab('Cell count') +
          ggtitle("Expected (red line) vs Simulated (blue dots) dynamics of total cell count")
}
