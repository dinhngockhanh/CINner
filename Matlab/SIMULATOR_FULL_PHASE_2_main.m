%==============================PHASE 2: COPY-NUMBER PROFILES OF A SAMPLE
function package_sample = SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
    global N_sample
    global size_CN_block_DNA N_chromosomes vec_CN_block_no
                                                                        tstart_total=tic;
%---------------------------------------------Input the clonal evolution
    T_current                       = package_clonal_evolution{1};
    N_clones                        = package_clonal_evolution{4};
    genotype_list_ploidy_chrom      = package_clonal_evolution{5};
    genotype_list_ploidy_block      = package_clonal_evolution{6};
    genotype_list_ploidy_allele     = package_clonal_evolution{7};
    evolution_traj_time             = package_clonal_evolution{13};
    evolution_traj_divisions        = package_clonal_evolution{14};
    evolution_traj_clonal_ID        = package_clonal_evolution{15};
    evolution_traj_population       = package_clonal_evolution{16};
%-------------------------------Find a random sample of final population
    final_clonal_ID                 = evolution_traj_clonal_ID{end};
    final_clonal_population         = evolution_traj_population{end};
    final_population                = [];
    for i=1:length(final_clonal_ID)
        clone                                           = final_clonal_ID(i);
        clonal_population                               = final_clonal_population(i);
        final_population(end+1:end+clonal_population)   = clone;
    end
    node_genotype_current           = datasample(final_population,N_sample,'Replace',false);
%---------------------------------Create CN object for the sampled cells
    sample_genotype                 = node_genotype_current;
%---Find the CN profiles for each clone found in the sample
    sample_genotype_unique          = unique(sample_genotype);
    sample_genotype_unique_profile  = cell(1,length(sample_genotype_unique));
    for i_clone=1:length(sample_genotype_unique)
%       Extract CN information for the clone from clonal evolution data
        clone_ID                    = sample_genotype_unique(i_clone);
        ploidy_chrom                = genotype_list_ploidy_chrom{clone_ID};
        ploidy_block                = genotype_list_ploidy_block{clone_ID};
        ploidy_allele               = genotype_list_ploidy_allele{clone_ID};
%       Build the CN profile in SIGNALS style for the clone
        vec_clone_chr               = [];
        vec_clone_start             = [];
        vec_clone_end               = [];
        vec_clone_copy              = [];
        vec_clone_state             = [];
        vec_clone_Min               = [];
        vec_clone_Maj               = [];
        for chrom=1:N_chromosomes
            chrom_block_count       = vec_CN_block_no(chrom);
            chrom_ploidy            = ploidy_chrom(chrom);
%           Find location information of each chromosome block
            vec_chr                 = chrom*ones(1,chrom_block_count);
            vec_start               = size_CN_block_DNA*(0:chrom_block_count-1)+1;
            vec_end                 = size_CN_block_DNA*(1:chrom_block_count);
%           Find CN counts for each allele of each chromosome block
            vec_Allele_1            = zeros(1,chrom_block_count);
            vec_Allele_2            = zeros(1,chrom_block_count);
            for strand=1:chrom_ploidy
                mat_allele                  = ploidy_allele{chrom,strand};
                for CN_row=1:size(mat_allele,1)
                    list_1                  = find(mat_allele(CN_row)==1);
                    vec_Allele_1(list_1)    = vec_Allele_1(list_1)+1;
                    list_2                  = find(mat_allele(CN_row)==2);
                    vec_Allele_2(list_2)    = vec_Allele_2(list_2)+1;
                end
            end
%           Find Major/Minor CN counts of each chromosome block
            if mean(vec_Allele_1)<=mean(vec_Allele_2)
                vec_Min             = vec_Allele_1;
                vec_Maj             = vec_Allele_2;
            else
                vec_Min             = vec_Allele_2;
                vec_Maj             = vec_Allele_1;
            end
%           Find total CN count of each chromosome block
            vec_copy                = vec_Min+vec_Maj;
            vec_state               = vec_copy;
%           Update the CN information of the clone
            vec_clone_chr           = [vec_clone_chr vec_chr];
            vec_clone_start         = [vec_clone_start vec_start];
            vec_clone_end           = [vec_clone_end vec_end];
            vec_clone_copy          = [vec_clone_copy vec_copy];
            vec_clone_state         = [vec_clone_state vec_state];
            vec_clone_Min           = [vec_clone_Min vec_Min];
            vec_clone_Maj           = [vec_clone_Maj vec_Maj];
        end
%       Store the CN profile for the clone
        genotype_unique_profile     = table(num2str(vec_clone_chr'),vec_clone_start',vec_clone_end',vec_clone_copy',vec_clone_state',vec_clone_Min',vec_clone_Maj',...
                                      'VariableNames',["chr","start","end","copy","state","Min","Maj"]);
        sample_genotype_unique_profile{i_clone} = genotype_unique_profile;
    end
%---Find the CN profiles for each cell in the sample
    sample_cell_ID                  = cell(1,N_sample);
    sample_clone_ID                 = sample_genotype;
    for i_cell=1:N_sample
        clone_ID                    = sample_genotype(i_cell);
        i_clone                     = find(sample_genotype_unique==clone_ID);
%       Find the CN profile for this cell
        cell_genotype_profile       = sample_genotype_unique_profile{i_clone};
%       Add column for cell ID
        cell_ID                     = ['Sample-Library-' num2str(i_cell) '-' num2str(i_cell)];

        sample_cell_ID{i_cell}      = cell_ID;

        vec_cell_id                 = cell(size(cell_genotype_profile,1),1);
        vec_cell_id(:)              = {cell_ID};
        cell_genotype_profile.cell_id   = vec_cell_id;
%       Update table of CN profiles for all cells in the sample
        if i_cell==1
            sample_genotype_profile = cell_genotype_profile;
        else
            sample_genotype_profile = [sample_genotype_profile;cell_genotype_profile];
        end
    end
%---------------------------------Output package of data from simulation
    package_sample{1}               = sample_genotype_profile;
    package_sample{2}               = sample_cell_ID;
    package_sample{3}               = sample_clone_ID;
                                                                        % fprintf('******************************************************************************************\n');
                                                                        % fprintf('* PHASE 2 (SAMPLE PHYLOGENY) - total runtime  =  %ds\n*\n',round(runtime_total));
                                                                        % fprintf('******************************************************************************************\n');
end
