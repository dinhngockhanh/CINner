%=================================================SIMULATE DRIVER EVENTS
function SIMULATOR_FULL_PHASE_1_drivers(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
    global genotype_list_driver_count genotype_list_driver_map
    global evolution_genotype_changes


    genotype_list_driver_count(genotype_daughter_1)         = genotype_list_driver_count(genotype_to_react);
    genotype_list_driver_map{genotype_daughter_1}           = genotype_list_driver_map{genotype_to_react};
    evolution_genotype_changes{genotype_daughter_1}{end+1}  = {};

    genotype_list_driver_count(genotype_daughter_2)         = genotype_list_driver_count(genotype_to_react);
    genotype_list_driver_map{genotype_daughter_2}           = genotype_list_driver_map{genotype_to_react};
    evolution_genotype_changes{genotype_daughter_2}{end+1}  = {};
end
