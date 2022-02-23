%==============================PHASE 3: COPY-NUMBER PROFILES OF A SAMPLE
function package_sample_phylogeny = SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample)
    global N_sample
%---------------------------------------------Input the clonal evolution
    T_current                                   = package_clonal_evolution{1};
    N_clones                                    = package_clonal_evolution{4};
    genotype_list_ploidy_chrom                  = package_clonal_evolution{5};
    genotype_list_ploidy_block                  = package_clonal_evolution{6};
    genotype_list_ploidy_allele                 = package_clonal_evolution{7};
    evolution_traj_time                         = package_clonal_evolution{13};
    evolution_traj_divisions                    = package_clonal_evolution{14};
    evolution_traj_clonal_ID                    = package_clonal_evolution{15};
    evolution_traj_population                   = package_clonal_evolution{16};
%-------------------------------------------------------Input the sample
    sample_cell_ID                              = package_sample{2};
    sample_clone_ID                             = package_sample{3};
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
%       Get time point
        time                                            = evolution_traj_time(i);
%       Get current clonal populations in total population
        eligible_clonal_ID                              = evolution_traj_clonal_ID{i+1};
        eligible_clonal_total_population                = evolution_traj_population{i+1};
%       Get current clonal populations in sample
        eligible_clonal_sample_population               = zeros(1,length(eligible_clonal_ID));
        for clone=1:length(eligible_clonal_ID)
            clone_ID                                    = eligible_clonal_ID(clone);
            eligible_clonal_sample_population(clone)    = length(find(node_genotype_current==clone_ID));
        end
%       Translate next clonal populations in total population as
%       thresholds that clonal populations in sample cannot exceed
        if i==1
            limit_clonal_total_population               = Inf(1,length(eligible_clonal_ID));
        else
            limit_clonal_total_population               = zeros(1,length(eligible_clonal_ID));
            eligible_clonal_ID_tmp                      = evolution_traj_clonal_ID{i};
            eligible_clonal_total_population_tmp        = evolution_traj_population{i};
            for clone=1:length(eligible_clonal_ID)
                clone_ID                                = eligible_clonal_ID(clone);
                loc_tmp                                 = find(eligible_clonal_ID_tmp==clone_ID);
                if ~isempty(loc_tmp)
                    limit_clonal_total_population(clone)= eligible_clonal_total_population_tmp(loc_tmp);
                end
            end
        end
%=======Sanity tests
        if sum(eligible_clonal_sample_population)~=length(node_genotype_current)
            fprintf('\nERROR: CLONAL POPULATIONS IN SAMPLE DO NOT ADD UP\n\n');
        elseif any(eligible_clonal_sample_population>eligible_clonal_total_population)
            fprintf('\nERROR: CLONAL POPULATIONS IN SAMPLE ARE LARGER THAN IN TOTAL CELL POPULATION\n\n');
            % eligible_clonal_ID
            % eligible_clonal_sample_population
            % eligible_clonal_total_population
            % evolution_traj_population{i+2}
            % disp('~~~~~~~~~~~~~~~~~~~~~')
            % mat_division_total_population
            % mat_division_sample
            % mat_division_sample_clone
            % disp('~~~~~~~~~~~~~~~~~~~~~')
            % limit_clonal_total_population
            % tmp_clonal_sample_population
            % disp('----------------------------------------------------------------------------')
        end
%=======Get list of divisions occurring in total population
%       Column 1:       number of divisions
%       Column 2:       genotype mother
%       Column 3:       genotype daughter 1
%       Column 4:       genotype daughter 2
        mat_division_total_population                   = evolution_traj_divisions{i};
        if isempty(mat_division_total_population)
            continue;
        end
%=======Reduce list of divisions to only clones present in sample
        mat_division_total_population_short             = [];
        for division=1:size(mat_division_total_population,1)
            clonal_ID_daughter_1                        = mat_division_total_population(division,3);
            clonal_population_daughter_1                = eligible_clonal_sample_population(find(eligible_clonal_ID==clonal_ID_daughter_1));
            clonal_ID_daughter_2                        = mat_division_total_population(division,4);
            clonal_population_daughter_2                = eligible_clonal_sample_population(find(eligible_clonal_ID==clonal_ID_daughter_2));
            if (clonal_population_daughter_1>0) || (clonal_population_daughter_2>0)
                mat_division_total_population_short     = [mat_division_total_population_short;mat_division_total_population(division,:)];
            end
        end
        mat_division_total_population                   = mat_division_total_population_short;
        if isempty(mat_division_total_population)
            continue;
        end
%=======One huge loop to make sure clonal populations in sample are correct
        logic_correct                                       = 0;
        while (logic_correct==0)
%-----------Simulate identities of all divisions occurring in sample
%           Row:            corresponding to mat_division_total_population
%           Column 1:       node indices undergoing division as daughter 1
%           Column 2:       node indices undergoing division as daughter 2
%           Column 3:       division indices for corresponding nodes on column 1
%           Column 4:       division indices for corresponding nodes on column 2
            mat_division_sample                             = cell(size(mat_division_total_population,1),4);
            for clone=1:length(eligible_clonal_ID)
%               For every clone found in the total population...
%               Find its clonal ID
                clonal_ID                                   = eligible_clonal_ID(clone);
%               Find its population in total population
                clonal_total_population                     = eligible_clonal_total_population(clone);
%               Find its population in sample's eligible nodes
                clonal_sample_population                    = eligible_clonal_sample_population(clone);
                if clonal_sample_population<=0
                    continue;
                end
%               Find all division roles that this clone plays in total population
%               Row 1:      division index (= row in mat_division_total_population)
%               Row 2:      daughter position (= 1 or 2)
%               Row 3:      Cell count for this division/position in total population
%               Row 4:      Node count for this division/position in sample ------> to be done in next sections
                mat_division_sample_clone                   = [];
                for daughter=1:2
                    vec_division_genotype_daughter          = mat_division_total_population(:,daughter+2);
                    vec_division                            = find(vec_division_genotype_daughter==clonal_ID);
                    mat_division_sample_clone_new           = [vec_division';daughter*ones(1,length(vec_division));mat_division_total_population(vec_division,1)'];
                    mat_division_sample_clone               = [mat_division_sample_clone mat_division_sample_clone_new];
                end
                if isempty(mat_division_sample_clone)
                    continue;
                end
%               Find total count of nodes of this clone to undergo divisions
%               of each type, i.e. row 4 in mat_division_sample_clone
                division_index_all                          = mat_division_sample_clone(1,:);
                count_nodes_each_max                        = mat_division_total_population(division_index_all,1)';
                freq                                        = sum(mat_division_sample_clone(3,:))/clonal_total_population;
%               Find total count of nodes to undergo divisions of all types
                count_nodes_all                             = binornd(clonal_sample_population,freq);
%               Divide total count of nodes among different division types
                count_nodes_each                            = mnrnd(count_nodes_all,mat_division_sample_clone(3,:)/sum(mat_division_sample_clone(3,:)));
                mat_division_sample_clone(4,:)              = count_nodes_each;
%               Check that node count in each position doesn't exceed limit in total population
                if any(count_nodes_each>count_nodes_each_max)
                    logic_correct                           = -1;
                    break;
                end
%               Jump to next clone if there is no division to perform
                if max(mat_division_sample_clone(4,:))==0
                    continue;
                end
%               Simulate which nodes undergo each division type
%               i.e. rows 1 & 2 in mat_division_sample
                eligible_nodes                              = node_list_current(find(node_genotype_current==clonal_ID));
                if length(eligible_nodes)==1
                    node_indices_all                        = eligible_nodes;
                else
                    node_indices_all                        = randsample(eligible_nodes,count_nodes_all);
                end
                for division_type=1:size(mat_division_sample_clone,2)
                    row                                     = mat_division_sample_clone(1,division_type);
                    col                                     = mat_division_sample_clone(2,division_type);
                    count                                   = mat_division_sample_clone(4,division_type);
                    if count>0
                        mat_division_sample{row,col}        = node_indices_all(1:count);
                    end
                    node_indices_all(1:count)               = [];
                end
%               Simulate the division indices for each division type
%               i.e. rows 3 & 4 in mat_division_sample
                for division_type=1:size(mat_division_sample_clone,2)
                    row                                     = mat_division_sample_clone(1,division_type);
                    col                                     = mat_division_sample_clone(2,division_type);
                    count                                   = mat_division_sample_clone(4,division_type);
                    if count>0
                        count_divisions_total               = mat_division_total_population(row,1);
                        mat_division_sample{row,col+2}      = transpose(randsample(count_divisions_total,count));
                    end
                end
            end
%           Redo the whole process if some node count in some position exceeded limit in total population
            if logic_correct==-1
                logic_correct                               = 0;
                continue;
            end
%-----------Update phylogeny tree according to the division identities
%           Save the current phylogeny in case new changes are wrong
            hclust_row_tmp                                  = hclust_row;
            hclust_nodes_tmp                                = hclust_nodes;
            hclust_merge_tmp                                = hclust_merge;
            hclust_height_tmp                               = hclust_height;

            phylogeny_origin_tmp                            = phylogeny_origin;
            phylogeny_elapsed_gens_tmp                      = phylogeny_elapsed_gens;
            phylogeny_elapsed_genotypes_tmp                 = phylogeny_elapsed_genotypes;
            phylogeny_genotype_tmp                          = phylogeny_genotype;
            phylogeny_birthtime_tmp                         = phylogeny_birthtime;
            phylogeny_deathtime_tmp                         = phylogeny_deathtime;

            node_list_current_tmp                           = node_list_current;
            node_genotype_current_tmp                       = node_genotype_current;
%           Update phylogeny according to every division
            for division_type=1:size(mat_division_total_population,1)
                genotype_mother                                     = mat_division_total_population(division_type,2);
%               Get list of nodes in positions of daughter 1 and daughter 2
                vec_nodes_daughter_1                                = mat_division_sample{division_type,1};
                vec_nodes_daughter_2                                = mat_division_sample{division_type,2};
                if isempty(vec_nodes_daughter_1) && isempty(vec_nodes_daughter_2)
                    continue;
                end
%               Get list of division indices of daughter 1 and daughter 2
                vec_div_indices_1                                   = mat_division_sample{division_type,3};
                vec_div_indices_2                                   = mat_division_sample{division_type,4};
                vec_div_indices_all                                 = unique([vec_div_indices_1 vec_div_indices_2]);
%               Perform each division
                for division=1:length(vec_div_indices_all)
                    div_index                                       = vec_div_indices_all(division);
                    loc_1                                           = find(vec_div_indices_1==div_index);
                    loc_2                                           = find(vec_div_indices_2==div_index);
                    if (~isempty(loc_1))&&(~isempty(loc_2))
%                       Nodes 1 and 2 are mergning...
                        node_1                                      = vec_nodes_daughter_1(loc_1);
                        node_2                                      = vec_nodes_daughter_2(loc_2);
                        node_mother                                 = min(node_list_current)-1;
%                       Update phylogeny in hclust style
                        hclust_row                                  = hclust_row+1;
                        hclust_nodes(node_mother)                   = hclust_row;
                        hclust_merge(hclust_row,:)                  = [hclust_nodes(node_1) hclust_nodes(node_2)];
                        hclust_height(hclust_row)                   = T_current-time;
%                       Update phylogeny in our style
                        phylogeny_origin(node_1)                    = node_mother;
                        phylogeny_origin(node_2)                    = node_mother;
                        phylogeny_elapsed_gens(node_mother)         = 1;
                        phylogeny_elapsed_genotypes{node_mother}    = [genotype_mother];
                        phylogeny_genotype(node_mother)             = genotype_mother;
                        phylogeny_birthtime(node_1)                 = time;
                        phylogeny_birthtime(node_2)                 = time;
                        phylogeny_deathtime(node_mother)            = time;
%                       Update phylogeny records in our style
                        pos_delete                                  = [find(node_list_current==node_1) find(node_list_current==node_2)];
                        node_list_current(pos_delete)               = [];
                        node_list_current                           = [node_mother node_list_current];
                        node_genotype_current(pos_delete)           = [];
                        node_genotype_current                       = [genotype_mother node_genotype_current];
                    else
%                       Either node 1 or node 2 has one more division...
                        if (~isempty(loc_1))
                            node_daughter                           = vec_nodes_daughter_1(loc_1);
                        else
                            node_daughter                           = vec_nodes_daughter_2(loc_2);
                        end
%                       Update phylogeny in our style
                        phylogeny_elapsed_gens(node_daughter)       = phylogeny_elapsed_gens(node_daughter)+1;
                        phylogeny_elapsed_genotypes{node_daughter}  = [genotype_mother phylogeny_elapsed_genotypes{node_daughter}];
%                       Update phylogeny records in our style
                        loc_daughter                                = find(node_list_current==node_daughter);
                        node_genotype_current(loc_daughter)         = genotype_mother;
                    end
                end
            end
%-----------Check if the clonal populations in sample satisfy conditions
%           Find clonal populations in sample after new divisions
            tmp_clonal_sample_population                    = zeros(1,length(eligible_clonal_ID));
            for clone=1:length(eligible_clonal_ID)
                clone_ID                                    = eligible_clonal_ID(clone);
                tmp_clonal_sample_population(clone)         = length(find(node_genotype_current==clone_ID));
            end
%           Redo this whole step if clonal populations violate thresholds
            if any(tmp_clonal_sample_population>limit_clonal_total_population)
                hclust_row                                  = hclust_row_tmp;
                hclust_nodes                                = hclust_nodes_tmp;
                hclust_merge                                = hclust_merge_tmp;
                hclust_height                               = hclust_height_tmp;

                phylogeny_origin                            = phylogeny_origin_tmp;
                phylogeny_elapsed_gens                      = phylogeny_elapsed_gens_tmp;
                phylogeny_elapsed_genotypes                 = phylogeny_elapsed_genotypes_tmp;
                phylogeny_genotype                          = phylogeny_genotype_tmp;
                phylogeny_birthtime                         = phylogeny_birthtime_tmp;
                phylogeny_deathtime                         = phylogeny_deathtime_tmp;

                node_list_current                           = node_list_current_tmp;
                node_genotype_current                       = node_genotype_current_tmp;
            else
                logic_correct                               = 1;
            end
        end
    end
%   Assign original cell to be born at the beginning of clonal evolution
    phylogeny_birthtime(1)                                  = evolution_traj_time(1);





    % hclust_row
    % hclust_nodes
    % hclust_labels
    %
    % hclust_merge
    % hclust_height
    %
    % phylogeny_origin
    % phylogeny_elapsed_gens
    % phylogeny_elapsed_genotypes
    % phylogeny_genotype
    % phylogeny_birthtime
    % phylogeny_deathtime
    %
    % node_genotype_current
    %
    % node_list_current







%-----------------------------------------Reorder the nodes for plotting
%---Find an order on all nodes of the phylogeny in our style
%   Find number of progeny of each node
    progeny_count                                   = zeros(1,2*N_sample-1);
    progeny_count(N_sample:end)                     = 1;
    for node=2*N_sample-1:-1:2
        mother_node                                 = phylogeny_origin(node);
        progeny_count(mother_node)                  = progeny_count(mother_node)+progeny_count(node);
    end
%   Reorder the sample phylogeny tree based on progeny counts
    phylogeny_order                                 = zeros(1,2*N_sample-1);
    phylogeny_order(1)                              = 1;
    for node=1:2*N_sample-1
        vec_daughter_nodes                          = find(phylogeny_origin==node);
        if length(vec_daughter_nodes)==2
            daughter_node_1                         = vec_daughter_nodes(1);
            progeny_count_1                         = progeny_count(daughter_node_1);
            daughter_node_2                         = vec_daughter_nodes(2);
            progeny_count_2                         = progeny_count(daughter_node_2);
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
    hclust_order_inverse                            = phylogeny_order(N_sample:end);
    hclust_order                                    = zeros(1,N_sample);
    for i_cell=1:N_sample
        loc                                         = hclust_order_inverse(i_cell);
        hclust_order(loc)                           = i_cell;
    end
%------------------------------------------------Create clustering table
%   Change clone ID from genotype index in simulation to [A,B,C,...]
    sample_clone_ID_numeric                         = sample_clone_ID;
    sample_clone_ID_unique_numeric                  = unique(sample_clone_ID_numeric);
    sample_clone_ID_unique_letters                  = {};
    for i=1:length(sample_clone_ID_unique_numeric)
        sample_clone_ID_unique_letters{i}           = char(i + 64);
    end
    sample_clone_ID_letters                         = cell(1,length(sample_clone_ID_numeric));
    for i_clone=1:length(sample_clone_ID_unique_numeric)
        clone_ID_numeric                            = sample_clone_ID_unique_numeric(i_clone);
        clone_ID_letters                            = sample_clone_ID_unique_letters{i_clone};
        vec_cell_ID                                 = find(sample_clone_ID_numeric==clone_ID_numeric);
        for i_cell=1:length(vec_cell_ID)
            cell_ID                                 = vec_cell_ID(i_cell);
            sample_clone_ID_letters{cell_ID}        = clone_ID_letters;
        end
    end
%   Create clustering table
    hclust_clustering                               = table(sample_cell_ID',sample_clone_ID_letters','VariableNames',["cell_id","clone_id"]);
%--------------------------------Create phylogeny object in hclust style
%   Create phylogeny object in hclust style
    phylogeny_hclust{1}                             = hclust_merge;
    phylogeny_hclust{2}                             = hclust_height;
    phylogeny_hclust{3}                             = hclust_labels;
    phylogeny_hclust{4}                             = hclust_order;
    phylogeny_hclust{5}                             = hclust_clustering;
%---------------------------------Output package of data from simulation
    package_sample_phylogeny{1}                     = phylogeny_hclust;
    package_sample_phylogeny{2}                     = phylogeny_origin;
    package_sample_phylogeny{3}                     = phylogeny_elapsed_gens;
    package_sample_phylogeny{4}                     = phylogeny_elapsed_genotypes;
    package_sample_phylogeny{5}                     = phylogeny_genotype;
    package_sample_phylogeny{6}                     = phylogeny_birthtime;
    package_sample_phylogeny{7}                     = phylogeny_deathtime;
    package_sample_phylogeny{8}                     = phylogeny_order;
end
