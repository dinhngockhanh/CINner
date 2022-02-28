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







print(clone_phylogeny_all_genotypes)

print(clone_hclust_nodes)

print(clone_hclust_merge)


    # print(clone_phylogeny_genotype)



}
