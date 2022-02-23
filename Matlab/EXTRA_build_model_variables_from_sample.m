%===============================PRODUCE MODEL VARIABLE FILES FROM SAMPLE
function EXTRA_build_model_variables_from_sample(model,folder_name,package_output)
    global N_chromosomes size_CN_block_DNA vec_CN_block_no vec_centromere_location
    global driver_library
%---------------------------------Input the clonal populations in sample
    package_clonal_evolution            = package_output{1};

    package_sample                      = package_output{2};

    genotype_list_ploidy_chrom          = package_clonal_evolution{5};
    genotype_list_ploidy_block          = package_clonal_evolution{6};
    genotype_list_ploidy_allele         = package_clonal_evolution{7};
    genotype_list_driver_count          = package_clonal_evolution{8};
    genotype_list_driver_map            = package_clonal_evolution{9};

    sample_clone_ID                     = package_sample{3};
%----------------------Output CN-driver profiles as model variable files
%---Find all unique clones in sample
    all_clones_ID                       = unique(sample_clone_ID);
    all_clones_population               = histc(sample_clone_ID,all_clones_ID);
    N_all_clones                        = length(all_clones_ID);
%---Create and output model variables - copy number state
    TABLE_CLONAL_COPY_NUMBER_PROFILES   = [];
    HEADER_CLONAL_COPY_NUMBER_PROFILES  = [];
%   Initialize table of CN profiles with chromosome and bin positions
    vec_chrom                           = [];
    vec_bin                             = [];
    for chrom=1:N_chromosomes
        bin_count                       = vec_CN_block_no(chrom);
        vec_chrom                       = [vec_chrom chrom*ones(1,bin_count)];
        vec_bin                         = [vec_bin [1:bin_count]];
    end
    TABLE_CLONAL_COPY_NUMBER_PROFILES   = num2cell([vec_chrom' vec_bin']);
    HEADER_CLONAL_COPY_NUMBER_PROFILES  = ["Chromosome","Bin"];
%   Update table of CN profiles with each clone found in sample
    for clone=1:N_all_clones
        clone_ID                    = all_clones_ID(clone);
%       Get the CN profile of this clone
        ploidy_chrom                = genotype_list_ploidy_chrom{clone_ID};
        ploidy_block                = genotype_list_ploidy_block{clone_ID};
        ploidy_allele               = genotype_list_ploidy_allele{clone_ID};
%       Find maximum strand count for any chromosome
        max_no_strands              = max(ploidy_chrom);
%       Translate CN profile of this clone into table
        TABLE_CLONE_CURRENT         = cell(size(TABLE_CLONAL_COPY_NUMBER_PROFILES,1),max_no_strands);
        N_row                       = 0;
        for chrom=1:N_chromosomes
            no_strands              = ploidy_chrom(chrom);
            bin_count               = vec_CN_block_no(chrom);
            for bin=1:bin_count
                N_row               = N_row+1;
                for strand=1:no_strands
                    unit_count      = ploidy_block{chrom,strand}(bin);
                    allele          = '';
                    for unit=1:unit_count
                        allele      = [allele char(ploidy_allele{chrom,strand}(unit,bin) + 64)];
                    end
if unit_count~=1
    fprintf(['Clone ' num2str(clone) ' - chrom ' num2str(chrom) ' - block ' num2str(bin) ' - strand ' num2str(strand) ': CN = ' num2str(unit_count) ' \n']);
end
                    if strcmp(allele,'')
                        allele      = 'NA';
                    end
                    TABLE_CLONE_CURRENT{N_row,strand}   = allele;
                end
                for strand=no_strands+1:max_no_strands
                    TABLE_CLONE_CURRENT{N_row,strand}   = 'NA';
                end
            end
        end
        TABLE_CLONAL_COPY_NUMBER_PROFILES   = [TABLE_CLONAL_COPY_NUMBER_PROFILES TABLE_CLONE_CURRENT];
%       Add in new headers
        HEADER_CLONE_CURRENT        = [];
        for strand=1:max_no_strands
            HEADER_CLONE_CURRENT    = [HEADER_CLONE_CURRENT string(['Clone_' num2str(clone) '_strand_' num2str(strand)])];
        end
        HEADER_CLONAL_COPY_NUMBER_PROFILES  = [HEADER_CLONAL_COPY_NUMBER_PROFILES HEADER_CLONE_CURRENT];
    end
%   Output model variables - copy number state
    TABLE_CLONAL_COPY_NUMBER_PROFILES   = cell2table(TABLE_CLONAL_COPY_NUMBER_PROFILES,'VariableNames',HEADER_CLONAL_COPY_NUMBER_PROFILES);
    filename                            = [folder_name model '-input-initial-cn-profiles.csv'];
    writetable(TABLE_CLONAL_COPY_NUMBER_PROFILES,filename);
%---Create and output model variables - other information
    vec_clone                           = [1:N_all_clones];
    vec_cell_count                      = all_clones_population;
    vec_drivers                         = [];
    for clone=1:N_all_clones
        DRIVERS                         = '';
        driver_map                      = genotype_list_driver_map{clone};
        for row=1:size(driver_map,1)
            driver_index                = driver_map(row,1);
            driver_ID                   = driver_library{driver_index,1}{1};
            driver_strand               = driver_map(row,3);
            driver_unit                 = driver_map(row,5);
            if ~isempty(DRIVERS)
                DRIVERS                 = [DRIVERS ';'];
            end
            DRIVERS                     = [DRIVERS driver_ID '_strand' num2str(driver_strand) '_unit' num2str(driver_unit)];
        end
        vec_drivers                     = [vec_drivers string(DRIVERS)];
    end
%   Output model variables - other information
    CLONAL_OTHERS                       = table(vec_clone',vec_cell_count',vec_drivers','VariableNames',["Clone","Cell_count","Drivers"]);
    filename                            = [folder_name model '-input-initial-others.csv'];
    writetable(CLONAL_OTHERS,filename);
end
