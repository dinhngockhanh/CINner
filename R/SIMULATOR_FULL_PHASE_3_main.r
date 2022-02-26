#==============================PHASE 3: COPY-NUMBER PROFILES OF A SAMPLE
SIMULATOR_FULL_PHASE_3_main <- function(package_clonal_evolution,package_sample) {
#---------------------------------------------Input the clonal evolution
    T_current                                               <- package_clonal_evolution[[1]]
    genotype_list_ploidy_chrom                              <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block                              <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele                             <- package_clonal_evolution[[7]]
    evolution_traj_time                                     <- package_clonal_evolution[[13]]
    evolution_traj_divisions                                <- package_clonal_evolution[[14]]
    evolution_traj_clonal_ID                                <- package_clonal_evolution[[15]]
    evolution_traj_population                               <- package_clonal_evolution[[16]]
#-------------------------------------------------------Input the sample
    sample_cell_ID                                          <- package_sample[[2]]
    sample_clone_ID                                         <- package_sample[[3]]
    sample_clone_ID_letters                                 <- package_sample[[4]]
    table_clone_ID_vs_letters                               <- package_sample[[5]]
    N_sample_clone                                          <- length(sample_clone_ID)
#-----------------------------------Initialize phylogeny in hclust style
#   Initialize information to build phylogeny in hclust style
    hclust_row                                              <- 0
    hclust_nodes                                            <- rep(0,1,2*N_sample-1)
    hclust_nodes[N_sample:(2*N_sample-1)]                   <- (-1:-N_sample)
    hclust_labels                                           <- sample_cell_ID
#   Initialize actual phylogeny in hclust style
    hclust_merge                                            <- matrix(0,nrow=N_sample-1,ncol=2)
    hclust_height                                           <- rep(0,1,N_sample-1)
#--------------------------------------Initialize phylogeny in our style
    phylogeny_origin                                        <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_gens                                  <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_genotypes                             <- vector("list",length=2*N_sample-1)
    phylogeny_genotype                                      <- rep(0,length=2*N_sample-1)
    phylogeny_birthtime                                     <- rep(0,length=2*N_sample-1)
    phylogeny_deathtime                                     <- rep(0,length=2*N_sample-1)
#   Initialize the current list of node genotypes
    node_genotype_current                                   <- sample_clone_ID
#   Initialize the current list of nodes in the sample phylogeny
    node_list_current                                       <- N_sample:(2*N_sample-1)
#   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens[node_list_current]               <- 1
    for (node in N_sample:(2*N_sample-1)) {
        phylogeny_elapsed_genotypes[[node]]                 <- node_genotype_current[node-N_sample+1]
    }
    phylogeny_genotype[node_list_current]                   <- node_genotype_current
    phylogeny_deathtime[node_list_current]                  <- T_current
#-----------------------------------Build the sample cell phylogeny tree
    for (i in seq(length(evolution_traj_divisions),1,-1)) {
    # for (i in seq(length(evolution_traj_divisions),length(evolution_traj_divisions),-1)) {
#       Get time point
        time                                                <- evolution_traj_time[i]
#       Get current clonal populations in total population
        eligible_clonal_ID                                  <- evolution_traj_clonal_ID[[i+1]]
        eligible_clonal_total_population                    <- evolution_traj_population[[i+1]]
#       Get current clonal populations in sample
        eligible_clonal_sample_population                   <- rep(0,length(eligible_clonal_ID))
        for (clone in 1:length(eligible_clonal_ID)) {
            clone_ID                                        <- eligible_clonal_ID[clone]
            eligible_clonal_sample_population[clone]        <- length(which(node_genotype_current==clone_ID))
        }
#       Translate next clonal populations in total population as
#       thresholds that clonal populations in sample cannot exceed
        if (i==1) {
            limit_clonal_total_population                   <- rep(Inf,length(eligible_clonal_ID))
            }
        else {
            limit_clonal_total_population                   <- rep(0,length(eligible_clonal_ID))
            eligible_clonal_ID_tmp                          <- evolution_traj_clonal_ID[[i]]
            eligible_clonal_total_population_tmp            <- evolution_traj_population[[i]]
            for (clone in 1:length(eligible_clonal_ID)) {
                clone_ID                                    <- eligible_clonal_ID[clone]
                loc_tmp                                     <- which(eligible_clonal_ID_tmp==clone_ID)
                if (length(loc_tmp)!=0) {
                    limit_clonal_total_population[clone]    <- eligible_clonal_total_population_tmp[loc_tmp]
                }
            }
        }
#=======Sanity tests
        if (sum(eligible_clonal_sample_population)!=length(node_genotype_current)) {
            fprintf('\nERROR: CLONAL POPULATIONS IN SAMPLE DO NOT ADD UP\n\n');
            }
        else { if (any(eligible_clonal_sample_population>eligible_clonal_total_population)) {
            cat('\nERROR: CLONAL POPULATIONS IN SAMPLE ARE LARGER THAN IN TOTAL CELL POPULATION\n\n')
            }
        }
#=======Get list of divisions occurring in total population
#       Column 1:       number of divisions
#       Column 2:       genotype mother
#       Column 3:       genotype daughter 1
#       Column 4:       genotype daughter 2
        mat_division_total_population                       <- evolution_traj_divisions[[i]]
        if (length(mat_division_total_population)==0) {
            next
        }
#=======One huge loop to make sure clonal populations in sample are correct
        logic_correct                                       <- 0
        while (logic_correct==0) {
#-----------Simulate identities of all divisions occurring in sample
#           Row:            corresponding to mat_division_total_population
#           Column 1:       node indices undergoing division as daughter 1
#           Column 2:       node indices undergoing division as daughter 2
#           Column 3:       division indices for corresponding nodes on column 1
#           Column 4:       division indices for corresponding nodes on column 2
#-----------Translation for R: linearized into vector:
#           Entries (1:4)   corresponds to row 1
#           Entries (5:8)   corresponds to row 2
#           ...
            mat_division_sample                             <- vector("list",length=(4*nrow(mat_division_total_population)))
            for (clone in 1:length(eligible_clonal_ID)) {
#               For every clone found in the total population...
#               Find its clonal ID
                clonal_ID                                   <- eligible_clonal_ID[clone]
#               Find its population in total population
                clonal_total_population                     <- eligible_clonal_total_population[clone]
#               Find its population in sample's eligible nodes
                clonal_sample_population                    <- eligible_clonal_sample_population[clone]
                if (clonal_sample_population<=0) {
                    next
                }
#               Find all division roles that this clone plays in total population
#               Row 1:      division index (= row in mat_division_total_population)
#               Row 2:      daughter position (= 1 or 2)
#               Row 3:      Cell count for this division/position in total population
#               Row 4:      Node count for this division/position in sample ------> to be done in next sections
                mat_division_sample_clone                   <- c()
                for (daughter in 1:2) {
                    vec_division_genotype_daughter          <- mat_division_total_population[,daughter+2]
                    vec_division                            <- which(vec_division_genotype_daughter==clonal_ID)
                    mat_division_sample_clone_new           <- rbind(matrix(vec_division,nrow=1),matrix(rep(daughter,length(vec_division)),nrow=1),matrix(mat_division_total_population[vec_division,1],nrow=1))
                    # mat_division_sample_clone_new           <- cbind(vec_division,rep(daughter,length(vec_division)),mat_division_total_population[vec_division,1])
                    mat_division_sample_clone               <- cbind(mat_division_sample_clone,mat_division_sample_clone_new)
                }
                if (length(mat_division_sample_clone)==0) {
                    next
                }
#               Find total count of nodes of this clone to undergo divisions
#               of each type, i.e. row 4 in mat_division_sample_clone
                division_index_all                          <- mat_division_sample_clone[1,]
                count_nodes_each_max                        <- mat_division_total_population[division_index_all,1]
                freq                                        <- sum(mat_division_sample_clone[3,])/clonal_total_population
#               Find total count of nodes to undergo divisions of all types
                count_nodes_all                             <- rbinom(n=1,size=clonal_sample_population,prob=freq)
#               Divide total count of nodes among different division types
                count_nodes_each                            <- rmultinom(n=1,size=count_nodes_all,mat_division_sample_clone[3,]/sum(mat_division_sample_clone[3,]))
                count_nodes_each                            <- matrix(count_nodes_each,nrow=1)
                mat_division_sample_clone                   <- rbind(mat_division_sample_clone,count_nodes_each)
#               Check that node count in each position doesn't exceed limit in total population
                if (any(count_nodes_each>count_nodes_each_max)) {
                    logic_correct                           <- -1
                    break
                }
#               Jump to next clone if there is no division to perform
                if (max(mat_division_sample_clone[4,])==0) {
                    next
                }
#               Simulate which nodes undergo each division type
#               i.e. columns 1 & 2 in mat_division_sample
                eligible_nodes                              <- node_list_current[which(node_genotype_current==clonal_ID)]
                if (length(eligible_nodes)==1) {
                    node_indices_all                        <- eligible_nodes
                    }
                else {
                    node_indices_all                        <- sample(x=eligible_nodes,size=count_nodes_all,replace=FALSE)
                }
                for (division_type in 1:ncol(mat_division_sample_clone)) {
                    row                                     <- mat_division_sample_clone[1,division_type]
                    col                                     <- mat_division_sample_clone[2,division_type]
                    count                                   <- mat_division_sample_clone[4,division_type]
                    if (count>0) {
                        ind                                 <- (row-1)*4+col
                        mat_division_sample[[ind]]          <- node_indices_all[1:count]
                        node_indices_all                    <- node_indices_all[-(1:count)]
                    }
                }
#               Simulate the division indices for each division type
#               i.e. columns 3 & 4 in mat_division_sample
                for (division_type in 1:ncol(mat_division_sample_clone)) {
                    row                                     <- mat_division_sample_clone[1,division_type]
                    col                                     <- mat_division_sample_clone[2,division_type]
                    count                                   <- mat_division_sample_clone[4,division_type]
                    if (count>0) {
                        col                                 <- col+2
                        ind                                 <- (row-1)*4+col
                        count_divisions_total               <- mat_division_total_population[row,1]
                        mat_division_sample[[ind]]          <- sample(x=count_divisions_total,size=count,replace=FALSE)
                    }
                }
            }
#           Redo the whole process if some node count in some position exceeded limit in total population
            if (logic_correct==-1) {
                logic_correct                               <- 0
                next
            }
#-----------Update phylogeny tree according to the division identities
#           Save the current phylogeny in case new changes are wrong
            hclust_row_tmp                                  <- hclust_row
            hclust_nodes_tmp                                <- hclust_nodes
            hclust_merge_tmp                                <- hclust_merge
            hclust_height_tmp                               <- hclust_height

            phylogeny_origin_tmp                            <- phylogeny_origin
            phylogeny_elapsed_gens_tmp                      <- phylogeny_elapsed_gens
            phylogeny_elapsed_genotypes_tmp                 <- phylogeny_elapsed_genotypes
            phylogeny_genotype_tmp                          <- phylogeny_genotype
            phylogeny_birthtime_tmp                         <- phylogeny_birthtime
            phylogeny_deathtime_tmp                         <- phylogeny_deathtime

            node_list_current_tmp                           <- node_list_current
            node_genotype_current_tmp                       <- node_genotype_current
#           Update phylogeny according to every division
            for (division_type in 1:nrow(mat_division_total_population)) {
                genotype_mother                                     <- mat_division_total_population[division_type,2]
#               Get list of nodes in positions of daughter 1 and daughter 2
                ind_1                                               <- (division_type-1)*4+1
                vec_nodes_daughter_1                                <- mat_division_sample[[ind_1]]
                ind_2                                               <- (division_type-1)*4+2
                vec_nodes_daughter_2                                <- mat_division_sample[[ind_2]]
                if ((length(vec_nodes_daughter_1)==0) && (length(vec_nodes_daughter_2)==0)) {
                    next
                }
#               Get list of division indices of daughter 1 and daughter 2
                ind_1                                               <- (division_type-1)*4+3
                vec_div_indices_1                                   <- mat_division_sample[[ind_1]]
                ind_2                                               <- (division_type-1)*4+4
                vec_div_indices_2                                   <- mat_division_sample[[ind_2]]
                vec_div_indices_all                                 <- unique(c(vec_div_indices_1,vec_div_indices_2))
#               Perform each division
                for (division in 1:length(vec_div_indices_all)) {
                    div_index                                       <- vec_div_indices_all[division]
                    loc_1                                           <- which(vec_div_indices_1==div_index)
                    loc_2                                           <- which(vec_div_indices_2==div_index)
                    if ((length(loc_1)!=0)&&(length(loc_2)!=0)) {
#                       Nodes 1 and 2 are mergning...
                        node_1                                      <- vec_nodes_daughter_1[loc_1]
                        node_2                                      <- vec_nodes_daughter_2[loc_2]
                        node_mother                                 <- min(node_list_current)-1
#                       Update phylogeny in hclust style
                        hclust_row                                  <- hclust_row+1
                        hclust_nodes[node_mother]                   <- hclust_row
                        hclust_merge[hclust_row,]                   <- c(hclust_nodes[node_1],hclust_nodes[node_2])
                        hclust_height[hclust_row]                   <- T_current-time
#                       Update phylogeny in our style
                        phylogeny_origin[node_1]                    <- node_mother
                        phylogeny_origin[node_2]                    <- node_mother
                        phylogeny_elapsed_gens[node_mother]         <- 1
                        phylogeny_elapsed_genotypes[[node_mother]]  <- c(genotype_mother)
                        phylogeny_genotype[node_mother]             <- genotype_mother
                        phylogeny_birthtime[node_1]                 <- time
                        phylogeny_birthtime[node_2]                 <- time
                        phylogeny_deathtime[node_mother]            <- time
#                       Update phylogeny records in our style
                        pos_delete                                  <- c(which(node_list_current==node_1),which(node_list_current==node_2))
                        node_list_current                           <- node_list_current[-pos_delete]
                        node_list_current                           <- c(node_mother,node_list_current)
                        node_genotype_current                       <- node_genotype_current[-pos_delete]
                        node_genotype_current                       <- c(genotype_mother,node_genotype_current)
                    }
                    else {
#                       Either node 1 or node 2 has one more division...
                        if (length(loc_1)!=0) {
                            node_daughter                           <- vec_nodes_daughter_1[loc_1]
                            }
                        else {
                            node_daughter                           <- vec_nodes_daughter_2[loc_2]
                        }
#                       Update phylogeny in our style
                        phylogeny_elapsed_gens[node_daughter]       <- phylogeny_elapsed_gens[node_daughter]+1
                        phylogeny_elapsed_genotypes[[node_daughter]]<- c(genotype_mother,phylogeny_elapsed_genotypes[[node_daughter]])
#                       Update phylogeny records in our style
                        loc_daughter                                <- which(node_list_current==node_daughter)
                        node_genotype_current[loc_daughter]         <- genotype_mother
                    }
                }
            }
#-----------Check if the clonal populations in sample satisfy conditions
#           Find clonal populations in sample after new divisions
            tmp_clonal_sample_population                            <- rep(0,length(eligible_clonal_ID))
            for (clone in 1:length(eligible_clonal_ID)) {
                clone_ID                                            <- eligible_clonal_ID[clone]
                tmp_clonal_sample_population[clone]                 <- length(which(node_genotype_current==clone_ID))
            }
#           Redo this whole step if clonal populations violate thresholds
            if (any(tmp_clonal_sample_population>limit_clonal_total_population)) {
                hclust_row                                  <- hclust_row_tmp
                hclust_nodes                                <- hclust_nodes_tmp
                hclust_merge                                <- hclust_merge_tmp
                hclust_height                               <- hclust_height_tmp

                phylogeny_origin                            <- phylogeny_origin_tmp
                phylogeny_elapsed_gens                      <- phylogeny_elapsed_gens_tmp
                phylogeny_elapsed_genotypes                 <- phylogeny_elapsed_genotypes_tmp
                phylogeny_genotype                          <- phylogeny_genotype_tmp
                phylogeny_birthtime                         <- phylogeny_birthtime_tmp
                phylogeny_deathtime                         <- phylogeny_deathtime_tmp

                node_list_current                           <- node_list_current_tmp
                node_genotype_current                       <- node_genotype_current_tmp
                }
            else {
                logic_correct                               <- 1
            }
        }
    }
#--------------------------------------------Complete the unmerged nodes
#-----------------------------i.e. there is more than one ancestral cell
#   Find all unmerged nodes
    list_unmerged_nodes                                     <- which(phylogeny_origin==0 & hclust_nodes!=0)
    list_unnecessary_nodes                                  <- which(phylogeny_origin==0 & hclust_nodes==0)
    N_unnecessary_nodes                                     <- length(list_unnecessary_nodes)
#---Complete the phylogeny in hclust style
    node_anchor                                             <- list_unmerged_nodes[1]
    hclust_node_anchor                                      <- hclust_nodes[node_anchor]
#   Merge all unmerged nodes together at first time point
    if (length(list_unmerged_nodes)>=2){
        for (i in 2:length(list_unmerged_nodes)){
            node                                            <- list_unmerged_nodes[i]
            hclust_row                                      <- hclust_row+1
            hclust_merge[hclust_row,]                       <- c(hclust_node_anchor,hclust_nodes[node])
            hclust_node_anchor                              <- hclust_node_anchor+1
            hclust_height[hclust_row]                       <- T_current
        }
    }
#---Complete the phylogeny in our style
#   Merge all unmerged nodes together at first time point
    phylogeny_birthtime[list_unmerged_nodes]                <- evolution_traj_time[1]
#   Delete unnecessary nodes
    phylogeny_origin                                        <- phylogeny_origin-N_unnecessary_nodes
    phylogeny_origin[list_unmerged_nodes]                   <- 0
    if (length(list_unnecessary_nodes)>0){
        phylogeny_origin                                    <- phylogeny_origin[-list_unnecessary_nodes]
        phylogeny_elapsed_gens                              <- phylogeny_elapsed_gens[-list_unnecessary_nodes]
        phylogeny_elapsed_genotypes                         <- phylogeny_elapsed_genotypes[-list_unnecessary_nodes]
        phylogeny_genotype                                  <- phylogeny_genotype[-list_unnecessary_nodes]
        phylogeny_birthtime                                 <- phylogeny_birthtime[-list_unnecessary_nodes]
        phylogeny_deathtime                                 <- phylogeny_deathtime[-list_unnecessary_nodes]
    }
#-----------------------------------------Reorder the nodes for plotting
    list_roots                                              <- list_unmerged_nodes-N_unnecessary_nodes
#---Find an order on all nodes of the phylogeny in our style
#   Find number of progeny of each node
    progeny_count                                           <- rep(0,length(phylogeny_origin))
    end                                                     <- length(progeny_count)
    progeny_count[(end-N_sample+1):end]                     <- 1
    for (node in length(progeny_count):1){
        mother_node                                         <- phylogeny_origin[node]
        if (mother_node>0){
            progeny_count[mother_node]                      <- progeny_count[mother_node]+progeny_count[node]
        }
    }
#   Reorder the sample phylogeny tree based on progeny counts
    phylogeny_order                                         <- rep(0,length(phylogeny_origin))
    phylogeny_order[list_roots]                             <- 1
    for (node in 0:length(progeny_count)){
        vec_daughter_nodes                                  <- which(phylogeny_origin==node)
        if (length(vec_daughter_nodes)==0){
            next
        }
        vec_progeny_counts                                  <- progeny_count[vec_daughter_nodes]
        tmp                                                 <- sort(vec_progeny_counts,index.return=TRUE)
        vec_progeny_counts                                  <- tmp$x
        vec_order                                           <- tmp$ix
        vec_daughter_nodes                                  <- vec_daughter_nodes[vec_order]
        for (i in 1:length(vec_daughter_nodes)){
            daughter_node                                   <- vec_daughter_nodes[i]
            if (i>1){
                progeny_count_extra                         <- sum(vec_progeny_counts[1:i-1])
            }else{
                progeny_count_extra                         <- 0
            }
            if (node==0){
                phylogeny_order[daughter_node]              <- phylogeny_order[daughter_node]+progeny_count_extra
            }
            else{
                phylogeny_order[daughter_node]              <- phylogeny_order[node]+progeny_count_extra
            }
        }
    }
#---Extract the order for phylogeny in hclust style
    hclust_order_inverse                                    <- phylogeny_order[(length(phylogeny_order)-N_sample+1):length(phylogeny_order)]
    hclust_order                                            <- rep(0,N_sample)
    for (i_cell in 1:N_sample){
        loc                                                 <- hclust_order_inverse[i_cell]
        hclust_order[loc]                                   <- i_cell
    }
#------------------------------------------------Create clustering table
    hclust_clustering                                       <- data.frame(sample_cell_ID,sample_clone_ID_letters)
    names(hclust_clustering)                                <- c('cell_id','clone_id')
#--------------------------------Create phylogeny object in hclust style
#   Create phylogeny object in hclust style
    phylogeny_hclust                                        <- list()
    phylogeny_hclust$merge                                  <- hclust_merge
    phylogeny_hclust$height                                 <- hclust_height
    phylogeny_hclust$order                                  <- hclust_order
    phylogeny_hclust$labels                                 <- sample_cell_ID
    class(phylogeny_hclust)                                 <- "hclust"






#-------------------------------
    clone_phylogeny_labels          <- table_clone_ID_vs_letters$Clone_ID_letter
    clone_phylogeny_ID              <- table_clone_ID_vs_letters$Clone_ID_number
#---Find cell MRCA node, merge time and genotypes for each clone
    clone_phylogeny_cell_MRCA       <- rep(0,length(clone_phylogeny_labels))
    clone_phylogeny_merge_time      <- rep(0,length(clone_phylogeny_labels))
    clone_phylogeny_genotypes       <- vector("list",length=length(clone_phylogeny_labels))

    vec_leaves_genotype                                             <- phylogeny_genotype
    vec_leaves_genotype[1:(length(phylogeny_genotype)-N_sample)]    <- -1
    for (clone in 1:length(clone_phylogeny_ID)){
        clone_ID                    <- clone_phylogeny_ID[clone]
#       Get node indices for cell leaves belonging in this clone
        vec_clone_leaves            <- which(vec_leaves_genotype==clone_ID)
#       Get list of all potential MRCA nodes for these cell leaves
        node                        <- vec_clone_leaves[1]
        node_mother                 <- phylogeny_origin[node]
        vec_potential_MRCA          <- c(node)
        while (node_mother>0){
            node                    <- node_mother
            node_mother             <- phylogeny_origin[node]
            vec_potential_MRCA      <- c(vec_potential_MRCA,node)
        }
#       Find MRCA node for all cell leaves belonging in this clone
        if (length(vec_clone_leaves)==1){
            node_MRCA               <- vec_clone_leaves
        }else{
            node_MRCA               <- vec_potential_MRCA[1]
            for (i in 2:length(vec_clone_leaves)){
                node                <- vec_clone_leaves[i]
                while (is.element(node,vec_potential_MRCA)==FALSE){
                    node            <- phylogeny_origin[node]
                }
                if (which(vec_potential_MRCA==node_MRCA)<which(vec_potential_MRCA==node)){
                    node_MRCA       <- node
                }
            }
        }
        clone_phylogeny_cell_MRCA[clone]    <- node_MRCA
        clone_phylogeny_merge_time[clone]   <- phylogeny_birthtime[node_MRCA]
#       Find all genotypes of this clone
        clone_genotypes             <- c()
        for (i in 1:length(vec_clone_leaves)){
            node                    <- vec_clone_leaves[i]

print('-----------------------------------')
print(node)
print(phylogeny_birthtime[node])
print(phylogeny_birthtime[node_MRCA])
            while (phylogeny_birthtime[node]>=phylogeny_birthtime[node_MRCA]){
            # while ((which(vec_potential_MRCA==node)<=(which(vec_potential_MRCA==node_MRCA))) & (node>0)){
# print('-----------------------------------')
# print(node)
# print(node_MRCA)
# print(vec_potential_MRCA)
# print(which(vec_potential_MRCA==node))
# print(which(vec_potential_MRCA==node_MRCA))


                clone_genotypes     <- unique(c(clone_genotypes,phylogeny_elapsed_genotypes[[node]]))
                node                <- phylogeny_origin[node]
            }
        }
        clone_phylogeny_genotypes[[clone]]  <- clone_genotypes
    }
print(clone_phylogeny_labels)
# clone_phylogeny_ID
# clone_phylogeny_cell_MRCA
print(clone_phylogeny_merge_time)
print(clone_phylogeny_genotypes)






#---------------------------------Create phylogeny object in phylo style
#   Create phylogeny object in phylo style
    phylogeny_phylo                                         <- ape::as.phylo(phylogeny_hclust,use.labels=TRUE)
#   Create object containing both phylo-style tree and clustering
    phylogeny_clustering_truth                              <- list()
    phylogeny_clustering_truth$tree                         <- phylogeny_phylo
    phylogeny_clustering_truth$clustering                   <- hclust_clustering
#---------------------------------Output package of data from simulation
    output                                                  <- list()
    output[[1]]                                             <- phylogeny_clustering_truth
    output[[2]]                                             <- phylogeny_origin
    output[[3]]                                             <- phylogeny_elapsed_gens
    output[[4]]                                             <- phylogeny_elapsed_genotypes
    output[[5]]                                             <- phylogeny_genotype
    output[[6]]                                             <- phylogeny_birthtime
    output[[7]]                                             <- phylogeny_deathtime
    output[[8]]                                             <- phylogeny_order
    return(output)
}
