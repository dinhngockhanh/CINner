#==================================================CREATE ONE SIMULATION
SIMULATOR_FULL_PROGRAM_one_simulation <- function(model,stage_final) {
    SIMULATOR_VARIABLES_for_simulation(model)
#------------------------------------------Simulate the clonal evolution
    flag_success                                    <- 0
    while (flag_success==0) {
        output                                      <- SIMULATOR_FULL_PHASE_1_main()
        flag_success                                <- output[[1]]
        package_clonal_evolution                    <- output[[2]]
    }
    if(stage_final==1){
        package_output                              <- list()
        package_output[[1]]                         <- package_clonal_evolution
        return(package_output)
    }

}
