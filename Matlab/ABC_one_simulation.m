%============================CREATE ONE SIMULATION WITH INPUT PARAMETERS
function package_output = ABC_one_simulation(model,parameter_set)
%---------------------------------------------Set up the model variables
    SIMULATOR_VARIABLES_for_fitting(model,parameter_set);
%------------------------------------------Simulate the clonal evolution
    flag_success                                    = 0;
    while flag_success==0
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   CLONAL EVOLUTION')
tic
        [flag_success,package_clonal_evolution]     = SIMULATOR_FULL_PHASE_1_main();
toc
    end
    package_output{1}                               = package_clonal_evolution;

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   SAMPLING')
tic
    package_sample                                  = SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution);
    package_output{2}                               = package_sample;
toc

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   SAMPLE PHYLOGENY')
tic
    package_sample_phylogeny                        = SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample);
    package_output{3}                               = package_sample_phylogeny;
toc

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   DONE WITH SIMULATION')

end
