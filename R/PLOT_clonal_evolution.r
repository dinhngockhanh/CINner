#=====================================PLOT CLONAL EVOLUTION AS FISH PLOT
PLOT_clonal_evolution <- function(package_simulation,vec_time_plot,unit){
#---------------------------------------------Input the clonal evolution
    package_clonal_evolution                        <- package_simulation[[1]]
    evolution_origin                                <- package_clonal_evolution[[11]]
#---------------------------------------------Input the clonal phylogeny
    package_sample_phylogeny                        <- package_simulation[[3]]
    package_clone_phylogeny                         <- package_sample_phylogeny[[4]]

    clone_phylogeny_labels                          <- package_clone_phylogeny[[1]]
    clone_phylogeny_genotype                        <- package_clone_phylogeny[[4]]
    clone_hclust_nodes                              <- package_clone_phylogeny[[7]]
    clone_hclust_merge                              <- package_clone_phylogeny[[8]]

    N_clones                                        <- length(clone_phylogeny_labels)





print('PLOT CLONAL EVOLUTION ........')
#----------------------------Build the genotype list for each clone node
    clone_phylogeny_all_genotypes                   <- vector('list',length=length(clone_phylogeny_genotype))
#   Initialize genotype lists for clone leaves
    for (clone in length(clone_phylogeny_genotype):(length(clone_phylogeny_genotype)-N_clones+1)){
        all_genotypes                               <- clone_phylogeny_genotype[clone]
        while (all_genotypes[1]!=0){
            ancestor_genotype                       <- evolution_origin[all_genotypes[1]]
            all_genotypes                           <- c(ancestor_genotype,all_genotypes)
        }
        clone_phylogeny_all_genotypes[[clone]]      <- all_genotypes
    }
#   Update genotype lists with information from clone hclust
    for (clone_hclust_mother_node in 1:nrow(clone_hclust_merge)){
#       Find mother clone node's index in phylogeny our style
        clone_phylogeny_mother_node                 <- which(is.element(clone_hclust_nodes,clone_hclust_mother_node))
#       Find daughter clone nodes' indices in phylogeny our style
        vec_clone_hclust_daughter_nodes             <- clone_hclust_merge[clone_hclust_mother_node,]
        vec_clone_phylogeny_daughter_nodes          <- which(is.element(clone_hclust_nodes,vec_clone_hclust_daughter_nodes))
#       Update genotype lists for the 3 clone nodes
        clone_phylogeny_daughter_node_1             <- vec_clone_phylogeny_daughter_nodes[[1]]
        clone_phylogeny_daughter_node_2             <- vec_clone_phylogeny_daughter_nodes[[2]]
        mother_genotypes                            <- intersect(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]],clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]])
        clone_phylogeny_all_genotypes[[clone_phylogeny_mother_node]]    <- mother_genotypes
        clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]]<- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]],mother_genotypes)
        clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]]<- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]],mother_genotypes)
    }






print(clone_phylogeny_all_genotypes)

print(clone_hclust_nodes)

print(clone_hclust_merge)


    # print(clone_phylogeny_genotype)



}
