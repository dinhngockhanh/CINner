%==============================================PHASE 2: SAMPLE PHYLOGENY
function [package_sample_phylogeny,runtime] = SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
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
%----------------------------------------Initialize the phylogeny record
    phylogeny_origin                = zeros(1,2*N_sample-1);
    phylogeny_elapsed_gens          = zeros(1,2*N_sample-1);
    phylogeny_elapsed_genotypes     = cell(1,2*N_sample-1);
    phylogeny_genotype              = zeros(1,2*N_sample-1);
    phylogeny_birthtime             = zeros(1,2*N_sample-1);
    phylogeny_deathtime             = zeros(1,2*N_sample-1);
%-------------------------------Find a random sample of final population
%   Initialize the current list of nodes in the sample phylogeny
    node_list_current               = [N_sample:2*N_sample-1];
%   Initialize the current list of genoytpes of the nodes
    final_clonal_ID                 = evolution_traj_clonal_ID{end};
    final_clonal_population         = evolution_traj_population{end};
    final_population                = [];
    for i=1:length(final_clonal_ID)
        clone                                           = final_clonal_ID(i);
        clonal_population                               = final_clonal_population(i);
        final_population(end+1:end+clonal_population)   = clone;
    end
    node_genotype_current           = datasample(final_population,N_sample,'Replace',false);
%   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens(node_list_current)   = 1;
    for node=N_sample:2*N_sample-1
        phylogeny_elapsed_genotypes{node}       = node_genotype_current(node-N_sample+1);
    end
    phylogeny_genotype(node_list_current)       = node_genotype_current;
    phylogeny_deathtime(node_list_current)      = T_current;



%---------------------------------Create CN object for the sampled cells
    sample_genotype                 = phylogeny_genotype(N_sample:2*N_sample-1);
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
    for i_cell=1:N_sample
        clone_ID                    = sample_genotype(i_cell);
        i_clone                     = find(sample_genotype_unique==clone_ID);
%       Find the CN profile for this cell
        cell_genotype_profile       = sample_genotype_unique_profile{i_clone};
%       Add column for cell ID
        cell_ID                     = ['Sample-Library-' num2str(i_cell) '-' num2str(i_cell)];
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
    package_sample_phylogeny{1}         = sample_genotype_profile;
return























%----------------------------------------Build the sample phylogeny tree
    for i=length(evolution_traj_divisions):-1:1
%       Get time point
        time                        = evolution_traj_time(i);
%       Report on progress
        % if (i==length(evolution_traj_divisions))||(floor(time/365)<floor(evolution_traj_time(i+1)/365))
        %     fprintf('Time = %d\n',round(time/365));
        % end
%       Get current total clonal population (after divisions)
        total_clonal_ID             = evolution_traj_clonal_ID{i+1};
        total_clonal_population     = evolution_traj_population{i+1};
%       Get current sample clonal population (after divisions)
        sample_clonal_population    = zeros(1,N_clones);
        for node=1:length(node_genotype_current)
            genotype                            = node_genotype_current(node);
            sample_clonal_population(genotype)  = sample_clonal_population(genotype)+1;
        end
%       Get list of eligible nodes of each genotype
        sample_eligible_nodes       = cell(1,N_clones);
        for node=1:length(node_genotype_current)
            genotype                            = node_genotype_current(node);
            sample_eligible_nodes{genotype}     = [sample_eligible_nodes{genotype} node_list_current(node)];
        end
%       Get list of divisions
        matrix_division             = evolution_traj_divisions{i};
%       For each type of divisions...
        for event_type=1:size(matrix_division,1)
%           Get number of divisions
            no_divisions        = matrix_division(event_type,1);
%           Get genotype of mother
            genotype_mother     = matrix_division(event_type,2);
%           Get genotype of 1st daughter
            genotype_daughter_1 = matrix_division(event_type,3);
            position_daughter_1 = find(total_clonal_ID==genotype_daughter_1);
%           Get genotype of 2nd daughter
            genotype_daughter_2 = matrix_division(event_type,4);
            position_daughter_2 = find(total_clonal_ID==genotype_daughter_2);
%           If daughter genotypes are not in current nodes, move on
            if (sample_clonal_population(genotype_daughter_1)<=0)&&(sample_clonal_population(genotype_daughter_2)<=0)
                continue
            end

%           For each specific division...
            for division=1:no_divisions
%               If these genotypes are not in current nodes, move on
                if (sample_clonal_population(genotype_daughter_1)<=0)&&(sample_clonal_population(genotype_daughter_2)<=0)
                    continue
                end
%               Choose the first daughter node
                logic_node_1    = rand<sample_clonal_population(genotype_daughter_1)/total_clonal_population(position_daughter_1);
                if logic_node_1==1
                    pos_node_1  = randi(sample_clonal_population(genotype_daughter_1));
                    node_1      = sample_eligible_nodes{genotype_daughter_1}(pos_node_1);
                    sample_eligible_nodes{genotype_daughter_1}(pos_node_1)  = [];
                    sample_clonal_population(genotype_daughter_1)           = sample_clonal_population(genotype_daughter_1)-1;
                    total_clonal_population(position_daughter_1)            = total_clonal_population(position_daughter_1)-1;
                else
                    node_1      = 0;
                    total_clonal_population(position_daughter_1)            = total_clonal_population(position_daughter_1)-1;
                end
%               Choose the second daughter node
                logic_node_2    = rand<sample_clonal_population(genotype_daughter_2)/total_clonal_population(position_daughter_2);
                if logic_node_2==1
                    pos_node_2  = randi(sample_clonal_population(genotype_daughter_2));
                    node_2      = sample_eligible_nodes{genotype_daughter_2}(pos_node_2);
                    sample_eligible_nodes{genotype_daughter_2}(pos_node_2)  = [];
                    sample_clonal_population(genotype_daughter_2)           = sample_clonal_population(genotype_daughter_2)-1;
                    total_clonal_population(position_daughter_2)            = total_clonal_population(position_daughter_2)-1;
                else
                    node_2      = 0;
                    total_clonal_population(position_daughter_2)            = total_clonal_population(position_daughter_2)-1;
                end
%               Update the nodes
                if (node_1==0)&&(node_2==0)
                    continue
                elseif (node_1>0)&&(node_2==0)
                    phylogeny_elapsed_gens(node_1)                          = phylogeny_elapsed_gens(node_1)+1;
                    phylogeny_elapsed_genotypes{node_1}                     = [genotype_mother phylogeny_elapsed_genotypes{node_1}];

                    node_genotype_current(find(node_list_current==node_1))  = genotype_mother;
                elseif (node_1==0)&&(node_2>0)
                    phylogeny_elapsed_gens(node_2)                          = phylogeny_elapsed_gens(node_2)+1;
                    phylogeny_elapsed_genotypes{node_2}                     = [genotype_mother phylogeny_elapsed_genotypes{node_2}];

                    node_genotype_current(find(node_list_current==node_2))  = genotype_mother;
                elseif (node_1>0)&&(node_2>0)
                    node_mother                                             = min(node_list_current)-1;

                    phylogeny_origin(node_1)                                = node_mother;
                    phylogeny_origin(node_2)                                = node_mother;
                    phylogeny_elapsed_gens(node_mother)                     = 1;
                    phylogeny_elapsed_genotypes{node_mother}                = [genotype_mother];
                    phylogeny_genotype(node_mother)                         = genotype_mother;
                    phylogeny_birthtime(node_1)                             = time;
                    phylogeny_birthtime(node_2)                             = time;
                    phylogeny_deathtime(node_mother)                        = time;

                    pos_delete                                              = [find(node_list_current==node_1) find(node_list_current==node_2)];
                    node_genotype_current(pos_delete)                       = [];
                    node_genotype_current                                   = [genotype_mother node_genotype_current];
                    node_list_current(pos_delete)                           = [];
                    node_list_current                                       = [node_mother node_list_current];
                end
            end
        end
    end
%---------------------------------Output package of data from simulation
    package_sample_phylogeny{1}         = phylogeny_origin;
    package_sample_phylogeny{2}         = phylogeny_elapsed_gens;
    package_sample_phylogeny{3}         = phylogeny_elapsed_genotypes;
    package_sample_phylogeny{4}         = phylogeny_genotype;
    package_sample_phylogeny{5}         = phylogeny_birthtime;
    package_sample_phylogeny{6}         = phylogeny_deathtime;
                                                                        runtime_total=toc(tstart_total);
    runtime{1}                          = runtime_total;
                                                                        % fprintf('******************************************************************************************\n');
                                                                        % fprintf('* PHASE 2 (SAMPLE PHYLOGENY) - total runtime  =  %ds\n*\n',round(runtime_total));
                                                                        % fprintf('******************************************************************************************\n');
end
