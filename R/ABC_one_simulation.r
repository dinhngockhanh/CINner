#============================CREATE ONE SIMULATION WITH INPUT PARAMETERS
ABC_one_simulation <- function(model,parameter_set) {
#---------------------------------------------Set up the model variables
    SIMULATOR_VARIABLES_for_fitting(model,parameter_set)
#------------------------------------------Simulate the clonal evolution
    flag_success                                    <- 0
    while (flag_success==0) {
        output                                      <- SIMULATOR_FULL_PHASE_1_main()
        flag_success                                <- output[[1]]
        package_clonal_evolution                    <- output[[2]]
    }


print(package_clonal_evolution)

}
