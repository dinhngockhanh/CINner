%==============================PHASE 3: COPY-NUMBER PROFILES OF A SAMPLE
function package_sample_phylogeny = SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample)
    global N_sample
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
%-------------------------------------------------------Input the sample
    sample_cell_ID                  = package_sample{2};
    sample_clone_ID                 = package_sample{3};










%-----------------------------------Initialize phylogeny in hclust style
%   Initialize information to build phylogeny in hclust style
    hclust_row                                  = 0;
    hclust_nodes                                = zeros(1,2*N_sample-1);
    hclust_nodes(N_sample:end)                  = -[1:N_sample];
    hclust_labels                               = sample_cell_ID;
%   Initialize actual phylogeny in hclust style
    hclust_merge                                = zeros(N_sample-1,2);
    hclust_height                               = zeros(1,N_sample-1);
%--------------------------------------Initialize phylogeny in our style
    phylogeny_origin                            = zeros(1,2*N_sample-1);
    phylogeny_elapsed_gens                      = zeros(1,2*N_sample-1);
    phylogeny_elapsed_genotypes                 = cell(1,2*N_sample-1);
    phylogeny_genotype                          = zeros(1,2*N_sample-1);
    phylogeny_birthtime                         = zeros(1,2*N_sample-1);
    phylogeny_deathtime                         = zeros(1,2*N_sample-1);
%   Initialize the current list of node genotypes
    node_genotype_current                       = sample_clone_ID;
%   Initialize the current list of nodes in the sample phylogeny
    node_list_current                           = [N_sample:2*N_sample-1];
%   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens(node_list_current)   = 1;
    for node=N_sample:2*N_sample-1
        phylogeny_elapsed_genotypes{node}       = node_genotype_current(node-N_sample+1);
    end
    phylogeny_genotype(node_list_current)       = node_genotype_current;
    phylogeny_deathtime(node_list_current)      = T_current;
%----------------------------------------Build the sample phylogeny tree
    for i=length(evolution_traj_divisions):-1:1
    % for i=length(evolution_traj_divisions):-1:length(evolution_traj_divisions)-10
%       Get time point
        time                                    = evolution_traj_time(i);
%       Report on progress
        total_clonal_ID                         = evolution_traj_clonal_ID{i+1};
        total_clonal_population                 = evolution_traj_population{i+1};
%       Get current sample clonal population (after divisions)
        sample_clonal_population                = zeros(1,N_clones);
        for node=1:length(node_genotype_current)
            genotype                            = node_genotype_current(node);
            sample_clonal_population(genotype)  = sample_clonal_population(genotype)+1;
        end
%       Get list of eligible nodes of each genotype
        sample_eligible_nodes                   = cell(1,N_clones);
        for node=1:length(node_genotype_current)
            genotype                            = node_genotype_current(node);
            sample_eligible_nodes{genotype}     = [sample_eligible_nodes{genotype} node_list_current(node)];
        end
%       Get list of divisions
        matrix_division                         = evolution_traj_divisions{i};
%       For each type of divisions...
        for event_type=1:size(matrix_division,1)
%           Get number of divisions
            no_divisions                        = matrix_division(event_type,1);
%           Get genotype of mother
            genotype_mother                     = matrix_division(event_type,2);
%           Get genotype of 1st daughter
            genotype_daughter_1                 = matrix_division(event_type,3);
            position_daughter_1                 = find(total_clonal_ID==genotype_daughter_1);
%           Get genotype of 2nd daughter
            genotype_daughter_2                 = matrix_division(event_type,4);
            position_daughter_2                 = find(total_clonal_ID==genotype_daughter_2);
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
                logic_node_1                                                = rand<sample_clonal_population(genotype_daughter_1)/total_clonal_population(position_daughter_1);
                if logic_node_1==1
                    pos_node_1                                              = randi(sample_clonal_population(genotype_daughter_1));
                    node_1                                                  = sample_eligible_nodes{genotype_daughter_1}(pos_node_1);
                    sample_eligible_nodes{genotype_daughter_1}(pos_node_1)  = [];
                    sample_clonal_population(genotype_daughter_1)           = sample_clonal_population(genotype_daughter_1)-1;
                    total_clonal_population(position_daughter_1)            = total_clonal_population(position_daughter_1)-1;
                else
                    node_1                                                  = 0;
                    total_clonal_population(position_daughter_1)            = total_clonal_population(position_daughter_1)-1;
                end
%               Choose the second daughter node
                logic_node_2                                                = rand<sample_clonal_population(genotype_daughter_2)/total_clonal_population(position_daughter_2);
                if logic_node_2==1
                    pos_node_2                                              = randi(sample_clonal_population(genotype_daughter_2));
                    node_2                                                  = sample_eligible_nodes{genotype_daughter_2}(pos_node_2);
                    sample_eligible_nodes{genotype_daughter_2}(pos_node_2)  = [];
                    sample_clonal_population(genotype_daughter_2)           = sample_clonal_population(genotype_daughter_2)-1;
                    total_clonal_population(position_daughter_2)            = total_clonal_population(position_daughter_2)-1;
                else
                    node_2                                                  = 0;
                    total_clonal_population(position_daughter_2)            = total_clonal_population(position_daughter_2)-1;
                end
%               Update the nodes
                if (node_1==0)&&(node_2==0)
%                   There is no merging....
                    continue
                elseif (node_1>0)&&(node_2==0)
%                   There is no merging but node 1 has one more division...
%                   Update phylogeny in our style
                    phylogeny_elapsed_gens(node_1)                          = phylogeny_elapsed_gens(node_1)+1;
                    phylogeny_elapsed_genotypes{node_1}                     = [genotype_mother phylogeny_elapsed_genotypes{node_1}];
%                   Update phylogeny records in our style
                    node_genotype_current(find(node_list_current==node_1))  = genotype_mother;
                elseif (node_1==0)&&(node_2>0)
%                   There is no merging but node 2 has one more division...
%                   Update phylogeny in our style
                    phylogeny_elapsed_gens(node_2)                          = phylogeny_elapsed_gens(node_2)+1;
                    phylogeny_elapsed_genotypes{node_2}                     = [genotype_mother phylogeny_elapsed_genotypes{node_2}];
%                   Update phylogeny records in our style
                    node_genotype_current(find(node_list_current==node_2))  = genotype_mother;
                elseif (node_1>0)&&(node_2>0)
%                   Nodes 1 and 2 are mergning...
                    node_mother                                             = min(node_list_current)-1;
%                   Update phylogeny in hclust style
                    hclust_row                                              = hclust_row+1;
                    hclust_nodes(node_mother)                               = hclust_row;
                    hclust_merge(hclust_row,:)                              = [hclust_nodes(node_1) hclust_nodes(node_2)];
                    hclust_height(hclust_row)                               = T_current-time;
%                   Update phylogeny in our style
                    phylogeny_origin(node_1)                                = node_mother;
                    phylogeny_origin(node_2)                                = node_mother;
                    phylogeny_elapsed_gens(node_mother)                     = 1;
                    phylogeny_elapsed_genotypes{node_mother}                = [genotype_mother];
                    phylogeny_genotype(node_mother)                         = genotype_mother;
                    phylogeny_birthtime(node_1)                             = time;
                    phylogeny_birthtime(node_2)                             = time;
                    phylogeny_deathtime(node_mother)                        = time;
%                   Update phylogeny records in our style
                    pos_delete                                              = [find(node_list_current==node_1) find(node_list_current==node_2)];
                    node_genotype_current(pos_delete)                       = [];
                    node_genotype_current                                   = [genotype_mother node_genotype_current];
                    node_list_current(pos_delete)                           = [];
                    node_list_current                                       = [node_mother node_list_current];
                end
            end
        end
    end
%   Assign original cell to be born at the beginning of clonal evolution
    phylogeny_birthtime(1)              = evolution_traj_time(1);




%-----------------------------------------Reorder the nodes for plotting
%---Find an order on all nodes of the phylogeny in our style
%   Find number of progeny of each node
    progeny_count                               = zeros(1,2*N_sample-1);
    progeny_count(N_sample:end)                 = 1;
    for node=2*N_sample-1:-1:2
        mother_node                             = phylogeny_origin(node);
        progeny_count(mother_node)              = progeny_count(mother_node)+progeny_count(node);
    end
%   Reorder the sample phylogeny tree based on progeny counts
    phylogeny_order                             = zeros(1,2*N_sample-1);
    phylogeny_order(1)                          = 1;
    for node=1:2*N_sample-1
        vec_daughter_nodes                      = find(phylogeny_origin==node);
        if length(vec_daughter_nodes)==2
            daughter_node_1                     = vec_daughter_nodes(1);
            progeny_count_1                     = progeny_count(daughter_node_1);
            daughter_node_2                     = vec_daughter_nodes(2);
            progeny_count_2                     = progeny_count(daughter_node_2);
            if progeny_count_1<progeny_count_2
                phylogeny_order(daughter_node_1)    = phylogeny_order(node);
                phylogeny_order(daughter_node_2)    = phylogeny_order(node)+progeny_count_1;
            else
                phylogeny_order(daughter_node_1)    = phylogeny_order(node)+progeny_count_2;
                phylogeny_order(daughter_node_2)    = phylogeny_order(node);
            end
        end
    end
%---Extract the order for phylogeny in hclust style
    hclust_order                                = phylogeny_order(N_sample:end);


%--------------------------------Create phylogeny object in hclust style
%   Create clustering table
    hclust_clustering                   = table(sample_cell_ID',sample_clone_ID','VariableNames',["cell_id","clone_id"]);
%   Create phylogeny object in hclust style
    hclust_object{1}                    = hclust_merge;
    hclust_object{2}                    = hclust_height;
    hclust_object{3}                    = hclust_labels;
    hclust_object{4}                    = hclust_order;
    hclust_object{5}                    = hclust_clustering;


%---------------------------------Output package of data from simulation
    package_sample_phylogeny{1}         = hclust_object;

    package_sample_phylogeny{2}         = phylogeny_origin;
    package_sample_phylogeny{3}         = phylogeny_elapsed_gens;
    package_sample_phylogeny{4}         = phylogeny_elapsed_genotypes;
    package_sample_phylogeny{5}         = phylogeny_genotype;
    package_sample_phylogeny{6}         = phylogeny_birthtime;
    package_sample_phylogeny{7}         = phylogeny_deathtime;
    package_sample_phylogeny{8}         = phylogeny_order;
end
