#===================OUTPUT: FOR COMPARISON BETWEEN MATLAB AND R VERSIONS
# function SIMULATOR_FULL_OUTPUT_for_comparison_between_MATLAB_and_R(folder_name,model,i_simulation,package_clonal_evolution,package_sample_phylogeny)
SIMULATOR_FULL_OUTPUT_for_comparison_between_MATLAB_and_R <- function(folder_name,model,i_simulation,package_clonal_evolution) {
#-----------------------------------------Input data from the simulation
    N_events                        <- package_clonal_evolution[[3]]
    N_clones                        <- package_clonal_evolution[[4]]
#-----------------------------------Output statistics for the comparison
    file_directory <- paste(folder_name,'/',model,'_sim_',i_simulation,'_statistics.csv',sep='')
    fileID         <- data.frame(N_clones,N_events)
    write.table(fileID,file_directory,row.names=FALSE,sep=',')
}
