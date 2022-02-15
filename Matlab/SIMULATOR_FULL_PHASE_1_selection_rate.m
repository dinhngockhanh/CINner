%==================================COMPUTE THE SELECTION RATE OF A CLONE
function clone_selection_rate = SIMULATOR_FULL_PHASE_1_selection_rate(driver_count,driver_map,ploidy_chrom,ploidy_block)
    global N_chromosomes vec_CN_block_no
    global bound_driver bound_ploidy
    global driver_library
%-------------------------Cell is not viable if losing whole chromosomes
%------------------------------------------or exceeding maximum CN count
    for chrom=1:N_chromosomes
        no_strands                  = ploidy_chrom(chrom);
        if no_strands<=0
            clone_selection_rate    = 0;
            return;
        end
        vec_CN                      = zeros(1,vec_CN_block_no(chrom));
        for strand=1:no_strands
            vec_CN                  = vec_CN+ploidy_block{chrom,strand};
        end
        if (max(vec_CN)==0) || (max(vec_CN)>bound_ploidy)
            clone_selection_rate    = 0;
            return;
        end
    end
%-------------------Cell is not viable if exceeding maximum driver count
    if driver_count>0
        driver_count_unique         = unique(driver_map(:,1));
        if length(driver_count_unique)>bound_driver
            clone_selection_rate    = 0;
            return;
        end
    end
%------------------------------Compute selection rates for viable clones
    driver_library_copy                             = driver_library;
    driver_library_copy.Copy_WT                     = zeros(size(driver_library,1),1);
    driver_library_copy.Copy_MUT                    = zeros(size(driver_library,1),1);
%   Find WT and MUT allele counts for each driver
    for i_driver=1:size(driver_library_copy,1)
        chrom                                       = driver_library_copy.Chromosome(i_driver);
        block                                       = driver_library_copy.Bin(i_driver);
        no_strands                                  = ploidy_chrom(chrom);
        driver_copy                                 = 0;
        for strand=1:no_strands
            driver_copy                             = driver_copy+ploidy_block{chrom,strand}(block);
        end
        driver_library_copy.Copy_WT(i_driver)       = driver_copy;
    end

    for i_driver=1:driver_count
        driver_ID                                   = driver_map(i_driver,1);
        driver_library_copy.Copy_MUT(driver_ID)     = driver_library_copy.Copy_MUT(driver_ID)+1;
        driver_library_copy.Copy_WT(driver_ID)      = driver_library_copy.Copy_WT(driver_ID)-1;
    end
%   Compute selection rate
    clone_selection_rate                            = prod(driver_library_copy.s_rate_WT.^driver_library_copy.Copy_WT)*prod(driver_library_copy.s_rate_MUT.^driver_library_copy.Copy_MUT);
end
