%============================CREATE ONE SIMULATION WITH INPUT PARAMETERS
function package_output = ABC_one_simulation(model,parameter_set)
%---------------------------------------------Set up the model variables
    SIMULATOR_VARIABLES_for_fitting(model,parameter_set);
%------------------------------------------Simulate the clonal evolution
    flag_success                                    = 0;
    while flag_success==0
        [flag_success,package_clonal_evolution]     = SIMULATOR_FULL_PHASE_1_main();
    end


    package_clonal_evolution


        % package_output{1}                           = package_clonal_evolution;









end
