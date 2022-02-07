%===================OUTPUT: FOR COMPARISON BETWEEN MATLAB AND R VERSIONS
function SIMULATOR_FULL_OUTPUT_for_comparison_between_MATLAB_and_R(folder_name,model,i_simulation,package_clonal_evolution)
%-----------------------------------------Input data from the simulation
    N_events_current                        = package_clonal_evolution{3};
    N_clones                                = package_clonal_evolution{4};
%-----------------------------------Output statistics for the comparison
    file_directory                          = [folder_name '/' model '_sim_' num2str(i_simulation) '_statistics.csv'];
    fileID                                  = fopen(file_directory,'w');
    fprintf(fileID,'N_clones,N_events\n%d,%d',N_clones,N_events_current);
end
