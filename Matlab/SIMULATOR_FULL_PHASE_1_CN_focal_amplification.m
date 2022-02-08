%===========================================SIMULATE FOCAL AMPLIFICATION
function SIMULATOR_FULL_PHASE_1_CN_focal_amplification(genotype_to_react,genotype_daughter)
    global genotype_list_ploidy_chrom genotype_list_ploidy_allele genotype_list_ploidy_block
    global genotype_list_driver_count genotype_list_driver_map genotype_list_DNA_length genotype_list_selection_rate
    global N_clones evolution_origin evolution_genotype_changes clonal_population_current evolution_traj_population

    global N_chromosomes size_CN_block_DNA vec_CN_block_no vec_centromere_location
    global growth_model bound_driver
    global driver_library

    global prob_CN_focal_amplified_length
%------------------------------------Find the new CN and driver profiles
%   Initialize the daughter's CN and driver profiles
    ploidy_chrom        = genotype_list_ploidy_chrom{genotype_daughter};
    ploidy_allele       = genotype_list_ploidy_allele{genotype_daughter};
    ploidy_block        = genotype_list_ploidy_block{genotype_daughter};
    driver_count        = genotype_list_driver_count(genotype_daughter);
    driver_map          = genotype_list_driver_map{genotype_daughter};
%   Find information about the focal amplification
    while 1
%       Choose the chromosome to be focally amplified
        chrom           = randi(N_chromosomes);
        no_strands      = ploidy_chrom(chrom);
        if no_strands<=0
            continue
        end
%       Choose the strand to be focally amplified
        strand          = randi(no_strands);
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
%       Choose the length of the focal amplification
        focal_length    = max_length+1;
        while focal_length>max_length
            focal_length= 1+geornd(prob_CN_focal_amplified_length);
        end
%       Choose the region to be focally amplified
        block_start     = (chrom_arm-1)*centromere + randi(max_length-focal_length+1);
        block_end       = block_start+focal_length-1;
        break;
    end
    % fprintf('\n=========================================================\nFOCAL AMPLIFICATION  ~~~  Chromosome %d, strand %d, block %d-%d\n',chrom,strand,block_start,block_end);
%   Find all drivers located on this region
    if (driver_count==0)||(isempty(intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end))))
        pos_drivers_to_amplify  = [];
    else
        pos_drivers_to_amplify  = intersect(intersect(find(driver_map(:,2)==chrom),find(driver_map(:,3)==strand)),intersect(find(driver_map(:,4)>=block_start),find(driver_map(:,4)<=block_end)));
    end
    N_drivers_to_amplify        = length(pos_drivers_to_amplify);
%   Update the chromosome strand allele identity
    for block=block_start:block_end
        block_CN                = ploidy_block{chrom,strand}(block);
        ploidy_allele{chrom,strand}(block_CN+1:2*block_CN,block)    = ploidy_allele{chrom,strand}(1:block_CN,block);
    end
%   Change the local CN on the amplified region
    ploidy_block{chrom,strand}(block_start:block_end)   = 2*ploidy_block{chrom,strand}(block_start:block_end);
%   Change the driver count
    driver_count                = driver_count+N_drivers_to_amplify;
%   Amplify the drivers
    if (N_drivers_to_amplify>0)
        driver_map(end+1:end+N_drivers_to_amplify,:)        = driver_map(pos_drivers_to_amplify,:);
        for driver=size(driver_map,1)+1-N_drivers_to_amplify:size(driver_map,1)
            block               = driver_map(driver,4);
            ploidy_block_driver = ploidy_block{chrom,strand}(block)/2;
            driver_map(driver,5)= driver_map(driver,5)+ploidy_block_driver;
        end
    end
%------------------------------------------------Output the new genotype
    genotype_list_ploidy_chrom{genotype_daughter}           = ploidy_chrom;
    genotype_list_ploidy_allele{genotype_daughter}          = ploidy_allele;
    genotype_list_ploidy_block{genotype_daughter}           = ploidy_block;
    genotype_list_driver_count(genotype_daughter)           = driver_count;
    genotype_list_driver_map{genotype_daughter}             = driver_map;
    evolution_genotype_changes{genotype_daughter}{end+1}    = {'focal-amplification',chrom,strand,block_start,block_end};
end
