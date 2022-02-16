#================================================CREATE MANY SIMULATIONS
SIMULATOR_FULL_PROGRAM_many_simulations <- function(model,stage_final,N_simulations) {
    for (i_simulation in 1:N_simulations){
        print('+++++++++++++++++++++++++++++++ SIMULATION',i_simulation)
        SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final)
    }
}
