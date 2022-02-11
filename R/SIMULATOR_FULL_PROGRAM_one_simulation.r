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
    if(stage_final>=1){
        package_output                              <- list()
        package_output[[1]]                         <- package_clonal_evolution
    }
#----------------------------------------------------Simulate the sample
    if(stage_final>=2){
        output                                      <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)

print('========================================================OUTPUT')
print(output)

        package_sample                              <- output[[1]]

print('========================================================package_sample')
print(package_sample)

        package_output[[2]]                         <- package_sample
    }
#-----------------------------------Simulate the phylogeny of the sample
    if(stage_final>=3){
        output                                      <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample)
        package_sample_phylogeny                    <- output[[1]]
        package_output[[3]]                         <- package_sample_phylogeny
    }
#------------------------------------------Output the simulation package
    return(package_output)
}
