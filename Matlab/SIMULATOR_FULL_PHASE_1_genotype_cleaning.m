%===================================DELETE NEW BUT UNNECESSARY GENOTYPES
function [genotype_to_react,genotype_daughter_1,genotype_daughter_2,position_to_react,position_daughter_1,position_daughter_2] = SIMULATOR_FULL_PHASE_1_genotype_cleaning(genotype_to_react,genotype_daughter_1,genotype_daughter_2,position_to_react,position_daughter_1,position_daughter_2)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length genotype_list_selection_rate genotype_list_prob_new_drivers
    global N_clones evolution_origin evolution_genotype_changes
    global clonal_population_current clonal_population_next clonal_ID_current
%-----------Find indicators of comparison between each pair of genotypes
    flag_mother_daughter1       = SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_to_react,genotype_daughter_1);
    flag_mother_daughter2       = SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_to_react,genotype_daughter_2);
    flag_daughter1_daughter2    = SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_daughter_1,genotype_daughter_2);
%--------------------------------------Delete new genotypes if necessary
    if (flag_mother_daughter1==0)&&(flag_mother_daughter2==0)&&(flag_daughter1_daughter2==0)
%       Case 1: all genotypes are different
        return;
    elseif (flag_mother_daughter1==0)&&(flag_mother_daughter2==1)&&(flag_daughter1_daughter2==0)
%       Case 2: combine genotype 2 into mother genotype
        N_clones                                                                = N_clones-1;

        evolution_origin(genotype_daughter_2)                                   = [];
        evolution_genotype_changes(genotype_daughter_2)                         = [];

        genotype_list_ploidy_chrom(genotype_daughter_2)                         = [];
        genotype_list_ploidy_allele(genotype_daughter_2)                        = [];
        genotype_list_ploidy_block(genotype_daughter_2)                         = [];
        genotype_list_driver_count(genotype_daughter_2)                         = [];
        genotype_list_driver_map(genotype_daughter_2)                           = [];
        genotype_list_DNA_length(genotype_daughter_2)                           = [];
        genotype_list_selection_rate(genotype_daughter_2)                       = [];
        genotype_list_prob_new_drivers(genotype_daughter_2)                     = [];

        clonal_ID_current(position_daughter_2)                                  = [];
        clonal_population_current(position_daughter_2)                          = [];
        clonal_population_next(position_daughter_2)                             = [];

        if genotype_daughter_1>genotype_daughter_2
            genotype_daughter_1                                                 = genotype_daughter_1-1;
        end
        genotype_daughter_2                                                     = genotype_to_react;

        if position_daughter_1>position_daughter_2
            position_daughter_1                                                 = position_daughter_1-1;
        end
        position_daughter_2                                                     = position_to_react;
    elseif (flag_mother_daughter1==1)&&(flag_mother_daughter2==0)&&(flag_daughter1_daughter2==0)
%       Case 3: combine genotype 1 into mother genotype
        N_clones                                                                = N_clones-1;

        evolution_origin(genotype_daughter_1)                                   = [];
        evolution_genotype_changes(genotype_daughter_1)                         = [];

        genotype_list_ploidy_chrom(genotype_daughter_1)                         = [];
        genotype_list_ploidy_allele(genotype_daughter_1)                        = [];
        genotype_list_ploidy_block(genotype_daughter_1)                         = [];
        genotype_list_driver_count(genotype_daughter_1)                         = [];
        genotype_list_driver_map(genotype_daughter_1)                           = [];
        genotype_list_DNA_length(genotype_daughter_1)                           = [];
        genotype_list_selection_rate(genotype_daughter_1)                       = [];
        genotype_list_prob_new_drivers(genotype_daughter_1)                     = [];

        clonal_ID_current(position_daughter_1)                                  = [];
        clonal_population_current(position_daughter_1)                          = [];
        clonal_population_next(position_daughter_1)                             = [];

        if genotype_daughter_2>genotype_daughter_1
            genotype_daughter_2                                                 = genotype_daughter_2-1;
        end
        genotype_daughter_1                                                     = genotype_to_react;

        if position_daughter_2>position_daughter_1
            position_daughter_2                                                 = position_daughter_2-1;
            clonal_ID_current(position_daughter_2)                              = clonal_ID_current(position_daughter_2)-1;
        end
        position_daughter_1                                                     = position_to_react;
    elseif (flag_mother_daughter1==0)&&(flag_mother_daughter2==0)&&(flag_daughter1_daughter2==1)
%       Case 4: combine genotype 2 into genotype 1
        N_clones                                                                = N_clones-1;

        evolution_origin(genotype_daughter_2)                                   = [];
        evolution_genotype_changes(genotype_daughter_2)                         = [];

        genotype_list_ploidy_chrom(genotype_daughter_2)                         = [];
        genotype_list_ploidy_allele(genotype_daughter_2)                        = [];
        genotype_list_ploidy_block(genotype_daughter_2)                         = [];
        genotype_list_driver_count(genotype_daughter_2)                         = [];
        genotype_list_driver_map(genotype_daughter_2)                           = [];
        genotype_list_DNA_length(genotype_daughter_2)                           = [];
        genotype_list_selection_rate(genotype_daughter_2)                       = [];
        genotype_list_prob_new_drivers(genotype_daughter_2)                     = [];

        clonal_ID_current(position_daughter_2)                                  = [];
        clonal_population_current(position_daughter_2)                          = [];
        clonal_population_next(position_daughter_2)                             = [];

        if genotype_daughter_1>genotype_daughter_2
            genotype_daughter_1                                                 = genotype_daughter_1-1;
        end
        genotype_daughter_2                                                     = genotype_daughter_1;

        if position_daughter_1>position_daughter_2
            position_daughter_1                                                 = position_daughter_1-1;
        end
        position_daughter_2                                                     = position_daughter_1;
    elseif (flag_mother_daughter1==1)&&(flag_mother_daughter2==1)&&(flag_daughter1_daughter2==1)
%       Case 5: combine both genotypes into mother genotype
        N_clones                                                                = N_clones-2;

        evolution_origin([genotype_daughter_1 genotype_daughter_2])             = [];
        evolution_genotype_changes([genotype_daughter_1 genotype_daughter_2])   = [];

        genotype_list_ploidy_chrom([genotype_daughter_1 genotype_daughter_2])   = [];
        genotype_list_ploidy_allele([genotype_daughter_1 genotype_daughter_2])  = [];
        genotype_list_ploidy_block([genotype_daughter_1 genotype_daughter_2])   = [];
        genotype_list_driver_count([genotype_daughter_1 genotype_daughter_2])   = [];
        genotype_list_driver_map([genotype_daughter_1 genotype_daughter_2])     = [];
        genotype_list_DNA_length([genotype_daughter_1 genotype_daughter_2])     = [];
        genotype_list_selection_rate([genotype_daughter_1 genotype_daughter_2]) = [];
        genotype_list_prob_new_drivers([genotype_daughter_1 genotype_daughter_2])= [];

        clonal_ID_current([position_daughter_1 position_daughter_2])            = [];
        clonal_population_current([position_daughter_1 position_daughter_2])    = [];
        clonal_population_next([position_daughter_1 position_daughter_2])       = [];

        genotype_daughter_1                                                     = genotype_to_react;
        genotype_daughter_2                                                     = genotype_to_react;

        position_daughter_1                                                     = position_to_react;
        position_daughter_2                                                     = position_to_react;
    end
end
