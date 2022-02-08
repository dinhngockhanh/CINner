%=================UPDATE DNA LENGTH AND SELECTION RATES OF NEW GENOTYPES
function SIMULATOR_FULL_PHASE_1_genotype_update(genotype_1,genotype_2)
    global rate_driver
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length
    global genotype_list_selection_rate genotype_list_prob_new_drivers
    global N_chromosomes size_CN_block_DNA
%---------------------------------Get the CN profile of the new genotype
    ploidy_chrom_1      = genotype_list_ploidy_chrom{genotype_1};
    ploidy_allele_1     = genotype_list_ploidy_allele{genotype_1};
    ploidy_block_1      = genotype_list_ploidy_block{genotype_1};
    driver_count_1      = genotype_list_driver_count(genotype_1);
    driver_map_1        = genotype_list_driver_map{genotype_1};

    ploidy_chrom_2      = genotype_list_ploidy_chrom{genotype_2};
    ploidy_allele_2     = genotype_list_ploidy_allele{genotype_2};
    ploidy_block_2      = genotype_list_ploidy_block{genotype_2};
    driver_count_2      = genotype_list_driver_count(genotype_2);
    driver_map_2        = genotype_list_driver_map{genotype_2};
%----------------------------Compute the DNA length for the new genotype
    DNA_length_1        = 0;
    for chrom=1:N_chromosomes
        for strand=1:ploidy_chrom_1(chrom)
            DNA_length_1= DNA_length_1+sum(ploidy_block_1{chrom,strand});
        end
    end
    DNA_length_1        = size_CN_block_DNA*DNA_length_1;
    DNA_length_2        = 0;
    for chrom=1:N_chromosomes
        for strand=1:ploidy_chrom_2(chrom)
            DNA_length_2= DNA_length_2+sum(ploidy_block_2{chrom,strand});
        end
    end
    DNA_length_2        = size_CN_block_DNA*DNA_length_2;
%-----------------------------Update the DNA length for the new genotype
    genotype_list_DNA_length{genotype_1}        = DNA_length_1;
    genotype_list_DNA_length{genotype_2}        = DNA_length_2;
%------------------------Compute the selection rate for the new genotype
    selection_rate_1    = SIMULATOR_FULL_PHASE_1_selection_rate(driver_count_1,driver_map_1,ploidy_chrom_1,ploidy_block_1);
    selection_rate_2    = SIMULATOR_FULL_PHASE_1_selection_rate(driver_count_2,driver_map_2,ploidy_chrom_2,ploidy_block_2);
%-----------------------------Update the DNA length for the new genotype
    genotype_list_selection_rate(genotype_1)    = selection_rate_1;
    genotype_list_selection_rate(genotype_2)    = selection_rate_2;
    genotype_list_prob_new_drivers(genotype_1)  = 1-poisspdf(0,rate_driver*DNA_length_1);
    genotype_list_prob_new_drivers(genotype_2)  = 1-poisspdf(0,rate_driver*DNA_length_2);
end
