%=================================COMPARE IF TWO GENOTYPES ARE IDENTICAL
function output = SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_1,genotype_2)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map
    global N_chromosomes
%------------------------------------------Initialize flag of comparison
    output          = 1;
%---------------------------------------Compare chromosome strand counts
    ploidy_chrom_1  = genotype_list_ploidy_chrom{genotype_1};
    ploidy_chrom_2  = genotype_list_ploidy_chrom{genotype_2};
    if ~isequal(ploidy_chrom_1,ploidy_chrom_2)
        output      = 0;
        return;
    end
%-----------------------------Compare block CN on each chromosome strand
    ploidy_block_1  = genotype_list_ploidy_block{genotype_1};
    ploidy_block_2  = genotype_list_ploidy_block{genotype_2};
    for chrom=1:N_chromosomes
        chrom_ploidy    = ploidy_chrom_1(chrom);
        for strand=1:chrom_ploidy
            if ~isequal(ploidy_block_1{chrom,strand},ploidy_block_2{chrom,strand})
                output  = 0;
                return;
            end
        end
    end
%--------------------------------------Compare chromosome strand alleles
    ploidy_allele_1 = genotype_list_ploidy_allele{genotype_1};
    ploidy_allele_2 = genotype_list_ploidy_allele{genotype_2};
    for chrom=1:N_chromosomes
        chrom_ploidy    = ploidy_chrom_1(chrom);
        for strand=1:chrom_ploidy
            if ~isequal(ploidy_allele_1{chrom,strand},ploidy_allele_2{chrom,strand})
                output  = 0;
                return;
            end
        end
    end
%--------------------------------------------------Compare driver counts
    driver_count_1  = genotype_list_driver_count(genotype_1);
    driver_count_2  = genotype_list_driver_count(genotype_2);
    if (driver_count_1~=driver_count_2)
        output  = 0;
        return;
    end
%----------------------------------------------------Compare driver maps
    driver_map_1    = genotype_list_driver_map{genotype_1};
    driver_map_2    = genotype_list_driver_map{genotype_2};
    if ~isequal(driver_map_1,driver_map_2)
        output  = 0;
        return;
    end
end
