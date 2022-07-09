#' @export
# ================================================CREATE MANY SIMULATIONS
ABC_many_simulations <- function(model, parameter_set, N_simulations) {
    for (i_simulation in 1:N_simulations) {
        print(i_simulation)
        ABC_one_simulation(model, parameter_set)
    }
}
