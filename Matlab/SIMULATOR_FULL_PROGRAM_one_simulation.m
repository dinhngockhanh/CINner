%==========================CREATE ONE SIMULATION WITH DEFAULT PARAMETERS
function package_output = SIMULATOR_FULL_PROGRAM_one_simulation(model,stage_final)
%---------------------------------------------Set up the model variables
    SIMULATOR_VARIABLES_for_simulation(model);
%------------------------------------------Simulate the clonal evolution
    flag_success                                    = 0;
    while flag_success==0
        [flag_success,package_clonal_evolution]     = SIMULATOR_FULL_PHASE_1_main();
    end
    if (stage_final>=1)
        package_output{1}                           = package_clonal_evolution;
    end
%----------------------------------------------------Simulate the sample
    if (stage_final>=2)
        package_sample                              = SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution);
    end
%-----------------------------------Simulate the phylogeny of the sample
    if (stage_final>=3)
        package_sample_phylogeny                    = SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample);
    end
end
