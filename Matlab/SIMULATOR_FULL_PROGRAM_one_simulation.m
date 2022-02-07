%==================================================CREATE ONE SIMULATION
function package_output = SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final)
    SIMULATOR_VARIABLES_for_simulation(model);
%------------------------------------------Simulate the clonal evolution
    flag_success                                    = 0;
    while flag_success==0
        [flag_success,package_clonal_evolution]     = SIMULATOR_FULL_PHASE_1_main();
    end
    if (stage_final==1)
        package_output{1}                           = package_clonal_evolution;
        return
    end

% %-------------------------------------Simulate the sample phylogeny tree
%     [package_sample_phylogeny,runtime_PHASE_2]                  = SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution);
% %------------------------------------------------Simulate the passengers
%     [package_passengers,runtime_PHASE_3]                        = SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample_phylogeny);
% %----------------------------------------------Simulate sampling effects
%     [package_SFS,runtime_PHASE_4]                               = SIMULATOR_FULL_PHASE_4_main(package_clonal_evolution,package_sample_phylogeny,package_passengers);
end
