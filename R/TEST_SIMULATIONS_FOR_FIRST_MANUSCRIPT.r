#=========================================SIMULATIONS FOR 1ST MANUSCRIPT
TEST_SIMULATIONS_FOR_FIRST_MANUSCRIPT <- function(model,stage_final,N_simulations) {
#-------------------------Create folder to store simulation output files
#     folder_name                 <- paste(model,'_R',sep='')
#     if (!file.exists(folder_name)) {
# #       If the folder with same name doesn't already exist, then create it
#         dir.create(folder_name)
#     }
#-----------------------------------------------------Create simulations
    package_output              <- list()
    for (i_simulation in 1:N_simulations) {
#-------Create one simulation
        cat('Simulation',i_simulation,'/',N_simulations,'\n')
#       Set up the variables for simulations
        SIMULATOR_VARIABLES_for_simulation(model)
#       Create one simulation
        output                              <- SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final)
        if (stage_final==1){
            package_clonal_evolution        <- output[[1]]
        }
        package_output[[i_simulation]]      <- package_clonal_evolution
#-------Statistics
#       Output statistics for comparison between MATLAB & R
        # SIMULATOR_FULL_OUTPUT_for_comparison_between_MATLAB_and_R(folder_name,model,i_simulation,package_clonal_evolution)
    }
    # return(package_output)
}
