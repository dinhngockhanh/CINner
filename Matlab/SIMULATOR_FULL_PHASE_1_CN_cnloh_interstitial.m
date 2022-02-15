%===============SIMULATE INTERSTITIAL COPY-NEUTRAL LOSS OF HETEROGENEITY
function SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial(genotype_to_react,genotype_daughter)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length genotype_list_selection_rate
    global N_clones evolution_origin evolution_genotype_changes clonal_population_current evolution_traj_population

    global N_chromosomes size_CN_block_DNA vec_CN_block_no vec_centromere_location
    global growth_model bound_driver
    global driver_library

    global prob_CN_cnloh_interstitial_length
%------------------------------------Find the new CN and driver profiles
%   Initialize the daughter's CN and driver profiles
    ploidy_chrom        = genotype_list_ploidy_chrom{genotype_daughter};
    ploidy_allele       = genotype_list_ploidy_allele{genotype_daughter};
    ploidy_block        = genotype_list_ploidy_block{genotype_daughter};
    driver_count        = genotype_list_driver_count(genotype_daughter);
    driver_map          = genotype_list_driver_map{genotype_daughter};
%   Find information about the interstitial CN-LOH
    while 1
%       Choose the chromosome to harbor the interstitial CN-LOH
        chrom           = randi(N_chromosomes);
        no_strands      = ploidy_chrom(chrom);
        if no_strands<=1
            continue
        end
%       Choose the strands to donate and receive the DNA
        vec_strands     = randsample(no_strands,2);
        strand_give     = vec_strands(1);
        strand_take     = vec_strands(2);
%       Find the chromosome's centromere location and length
        centromere      = vec_centromere_location(chrom);
        chrom_length    = vec_CN_block_no(chrom);
%       Choose the chromosome arm to be focally amplified
        chrom_arm       = randi(2);
        if chrom_arm==1
            max_length  = centromere;
        elseif chrom_arm==2
            max_length  = chrom_length-centromere;
        end
%       Choose the length of the interstitial CN-LOH
        cnloh_length    = max_length+1;
        while cnloh_length>max_length
            cnloh_length= 1+geornd(prob_CN_cnloh_interstitial_length);
        end
%       Choose the region to harbor the interstitial CN-LOH
        block_start     = (chrom_arm-1)*centromere + randi(max_length-cnloh_length+1);
        block_end       = block_start+cnloh_length-1;
        break;
    end
    % fprintf('\n=========================================================\nINTERSTITIAL CN-LOH  ~~~  Chromosome %d, strand %d to strand %d, block %d-%d\n',chrom,strand_give,strand_take,block_start,block_end);
%   Find all drivers to lose in the strand to receive the DNA
    if (driver_count==0)||(isempty(intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end))))
        pos_drivers_to_delete   = [];
    else
        pos_drivers_to_delete   = intersect(intersect(find(driver_map(:,2)==chrom),find(driver_map(:,3)==strand_take)),intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end)));
    end
    N_drivers_to_delete         = length(pos_drivers_to_delete);
%   Find all drivers to gain in the strand to donate the DNA
    if (driver_count==0)||(isempty(intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end))))
        pos_drivers_to_gain     = [];
    else
        pos_drivers_to_gain     = intersect(intersect(find(driver_map(:,2)==chrom),find(driver_map(:,3)==strand_give)),intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end)));
    end
    N_drivers_to_gain           = length(pos_drivers_to_gain);
%   Update the chromosome strand allele identity of the strand to receive the DNA
    ploidy_allele_strand_give   = ploidy_allele{chrom,strand_give};
    ploidy_allele{chrom,strand_take}(1:end,block_start:block_end)                               = 0;
    ploidy_allele{chrom,strand_take}(1:size(ploidy_allele_strand_give,1),block_start:block_end) = ploidy_allele_strand_give(1:end,block_start:block_end);
%   Change the local CN on the strand to receive the DNA
    ploidy_block{chrom,strand_take}(block_start:block_end)                                      = ploidy_block{chrom,strand_give}(block_start:block_end);
%   Copy the drivers gained during interstitial CN-LOH
    if (N_drivers_to_gain>0)
        driver_map_new                      = driver_map(pos_drivers_to_gain,:);
        driver_map_new(:,3)                 = strand_take;
        driver_map                          = [driver_map;driver_map_new];
    end
%   Remove the drivers lost during interstitial CN-LOH
    if (N_drivers_to_delete>0)
        driver_map(pos_drivers_to_delete,:) = [];
    end
%   Change the driver count
    driver_unique                           = unique(driver_map(:,1));
    driver_unique                           = driver_unique(driver_unique~=0);
    driver_count                            = length(driver_unique);
%------------------------------------------------Output the new genotype
    genotype_list_ploidy_chrom{genotype_daughter}           = ploidy_chrom;
    genotype_list_ploidy_allele{genotype_daughter}          = ploidy_allele;
    genotype_list_ploidy_block{genotype_daughter}           = ploidy_block;
    genotype_list_driver_count(genotype_daughter)           = driver_count;
    genotype_list_driver_map{genotype_daughter}             = driver_map;
    evolution_genotype_changes{genotype_daughter}{end+1}    = {'cnloh-interstitial',chrom,strand_give,strand_take,block_start,block_end};
end
