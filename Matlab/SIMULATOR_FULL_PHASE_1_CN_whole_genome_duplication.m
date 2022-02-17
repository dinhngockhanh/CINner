%======================================SIMULATE WHOLE GENOME DUPLICATION
function SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length genotype_list_selection_rate
    global N_clones evolution_origin evolution_genotype_changes clonal_population_current evolution_traj_population

    global N_chromosomes size_CN_block_DNA vec_CN_block_no
    global growth_model bound_driver
    global driver_library
%------------------------------------Find the new CN and driver profiles
%   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1      = genotype_list_ploidy_chrom{genotype_daughter_1};
    ploidy_allele_1     = genotype_list_ploidy_allele{genotype_daughter_1};
    ploidy_block_1      = genotype_list_ploidy_block{genotype_daughter_1};
    driver_count_1      = genotype_list_driver_count(genotype_daughter_1);
    driver_map_1        = genotype_list_driver_map{genotype_daughter_1};

    ploidy_chrom_2      = genotype_list_ploidy_chrom{genotype_daughter_2};
    ploidy_allele_2     = genotype_list_ploidy_allele{genotype_daughter_2};
    ploidy_block_2      = genotype_list_ploidy_block{genotype_daughter_2};
    driver_count_2      = genotype_list_driver_count(genotype_daughter_2);
    driver_map_2        = genotype_list_driver_map{genotype_daughter_2};
%   Change the chromosome ploidy of daughter cells
    ploidy_chrom_1      = 2*ploidy_chrom_1;
    ploidy_chrom_2      = 2*ploidy_chrom_2;
%   Update the chromosome strand allele identities of daughter cells
    for chrom=1:N_chromosomes
        chrom_ploidy                                        = ploidy_chrom_1(chrom);
        for strand=1:chrom_ploidy/2
            ploidy_allele_1{chrom,chrom_ploidy/2+strand}    = ploidy_allele_1{chrom,strand};
        end
        chrom_ploidy                                        = ploidy_chrom_2(chrom);
        for strand=1:chrom_ploidy/2
            ploidy_allele_2{chrom,chrom_ploidy/2+strand}    = ploidy_allele_2{chrom,strand};
        end
    end
%   Multiply the chromosome strands in each daughter cell
    for chrom=1:N_chromosomes
        chrom_ploidy                                        = ploidy_chrom_1(chrom);
        for strand=1:chrom_ploidy/2
            ploidy_block_1{chrom,chrom_ploidy/2+strand}     = ploidy_block_1{chrom,strand};
        end
        chrom_ploidy                                        = ploidy_chrom_2(chrom);
        for strand=1:chrom_ploidy/2
            ploidy_block_2{chrom,chrom_ploidy/2+strand}     = ploidy_block_2{chrom,strand};
        end
    end
%   Multiply the drivers in each daughter cell
    if driver_count_1>0
        driver_map_new_1                                    = driver_map_1;
        for driver=1:size(driver_map_new_1,1)
            chrom                                           = driver_map_1(driver,2);
            chrom_ploidy                                    = ploidy_chrom_1(chrom);
            driver_map_new_1(driver,3)                      = driver_map_new_1(driver,3)+chrom_ploidy/2;
        end
        driver_map_1                                        = [driver_map_1;driver_map_new_1];
    end
    if driver_count_2>0
        driver_map_new_2                                    = driver_map_2;
        for driver=1:size(driver_map_new_2,1)
            chrom                                           = driver_map_2(driver,2);
            chrom_ploidy                                    = ploidy_chrom_2(chrom);
            driver_map_new_2(driver,3)                      = driver_map_new_2(driver,3)+chrom_ploidy/2;
        end
        driver_map_2                                        = [driver_map_2;driver_map_new_2];
    end
%   Change the driver count in each daughter cell
    driver_unique_1                                         = unique(driver_map_1(:,1));
    driver_unique_1                                         = driver_unique_1(driver_unique_1~=0);
    driver_count_1                                          = length(driver_unique_1);

    driver_unique_2                                         = unique(driver_map_2(:,1));
    driver_unique_2                                         = driver_unique_2(driver_unique_2~=0);
    driver_count_2                                          = length(driver_unique_2);
%-----------------------------------------------Output the new genotypes
    genotype_list_ploidy_chrom{genotype_daughter_1}         = ploidy_chrom_1;
    genotype_list_ploidy_allele{genotype_daughter_1}        = ploidy_allele_1;
    genotype_list_ploidy_block{genotype_daughter_1}         = ploidy_block_1;
    genotype_list_driver_count(genotype_daughter_1)         = driver_count_1;
    genotype_list_driver_map{genotype_daughter_1}           = driver_map_1;
    evolution_genotype_changes{genotype_daughter_1}{end+1}  = {'whole-genome-duplication'};

    genotype_list_ploidy_chrom{genotype_daughter_2}         = ploidy_chrom_2;
    genotype_list_ploidy_allele{genotype_daughter_2}        = ploidy_allele_2;
    genotype_list_ploidy_block{genotype_daughter_2}         = ploidy_block_2;
    genotype_list_driver_count(genotype_daughter_2)         = driver_count_2;
    genotype_list_driver_map{genotype_daughter_2}           = driver_map_2;
    evolution_genotype_changes{genotype_daughter_2}{end+1}  = {'whole-genome-duplication'};
end
