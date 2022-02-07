%===========================INITIATE NEW GENOTYPES AFTER DRIVER/CN EVENT
function [genotype_daughter_1,genotype_daughter_2,position_daughter_1,position_daughter_2] = SIMULATOR_FULL_PHASE_1_genotype_initiation(genotype_to_react)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length genotype_list_selection_rate genotype_list_prob_new_drivers
    global N_clones evolution_origin evolution_genotype_changes clonal_population_current clonal_population_next clonal_ID_current
%----------------------Create a new genotype for the first daughter cell
    N_clones                                        = N_clones+1;
    genotype_daughter_1                             = N_clones;

    position_daughter_1                             = length(clonal_population_current)+1;
    clonal_ID_current(position_daughter_1)          = N_clones;
    clonal_population_current(position_daughter_1)  = 0;
    clonal_population_next(position_daughter_1)     = 0;

    evolution_origin(N_clones)                      = genotype_to_react;
    evolution_genotype_changes{N_clones}            = {};

    genotype_list_ploidy_chrom{N_clones}            = genotype_list_ploidy_chrom{genotype_to_react};
    genotype_list_ploidy_allele{N_clones}           = genotype_list_ploidy_allele{genotype_to_react};
    genotype_list_ploidy_block{N_clones}            = genotype_list_ploidy_block{genotype_to_react};
    genotype_list_driver_count(N_clones)            = genotype_list_driver_count(genotype_to_react);
    genotype_list_driver_map{N_clones}              = genotype_list_driver_map{genotype_to_react};
    genotype_list_DNA_length{N_clones}              = 0;
    genotype_list_selection_rate(N_clones)          = 0;
    genotype_list_prob_new_drivers(N_clones)        = 0;
%---------------------Create a new genotype for the second daughter cell
    N_clones                                        = N_clones+1;
    genotype_daughter_2                             = N_clones;

    position_daughter_2                             = length(clonal_population_current)+1;
    clonal_ID_current(position_daughter_2)          = N_clones;
    clonal_population_current(position_daughter_2)  = 0;
    clonal_population_next(position_daughter_2)     = 0;

    evolution_origin(N_clones)                      = genotype_to_react;
    evolution_genotype_changes{N_clones}            = {};

    genotype_list_ploidy_chrom{N_clones}            = genotype_list_ploidy_chrom{genotype_to_react};
    genotype_list_ploidy_allele{N_clones}           = genotype_list_ploidy_allele{genotype_to_react};
    genotype_list_ploidy_block{N_clones}            = genotype_list_ploidy_block{genotype_to_react};
    genotype_list_driver_count(N_clones)            = genotype_list_driver_count(genotype_to_react);
    genotype_list_driver_map{N_clones}              = genotype_list_driver_map{genotype_to_react};
    genotype_list_DNA_length{N_clones}              = 0;
    genotype_list_selection_rate(N_clones)          = 0;
    genotype_list_prob_new_drivers(N_clones)        = 0;
end
