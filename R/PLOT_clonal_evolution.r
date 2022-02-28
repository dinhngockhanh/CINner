#=====================================PLOT CLONAL EVOLUTION AS FISH PLOT
PLOT_clonal_evolution <- function(package_simulation,vec_time_plot,unit){
#---------------------------------------------Input the clonal evolution
    package_clonal_evolution                        <- package_simulation[[1]]
    evolution_origin                                <- package_clonal_evolution[[11]]
#---------------------------------------------Input the clonal phylogeny
    package_sample_phylogeny                        <- package_simulation[[3]]
    clone_phylogeny_labels                          <- package_sample_phylogeny[[10]]
    clone_phylogeny_genotype                        <- package_sample_phylogeny[[13]]

    N_clones                                        <- length(clone_phylogeny_labels)





print('PLOT CLONAL EVOLUTION ........')
#----------------------------Build the genotype list for each clone node
    clone_phylogeny_all_genotypes                   <- vector('list',length=length(clone_phylogeny_genotype))
#   Initialize genotype lists for clone leaves
    for (clone in length(clone_phylogeny_genotype):(length(clone_phylogeny_genotype)-N_clones+1)){
        clone_phylogeny_all_genotypes[[clone]]      <- clone_phylogeny_genotype[clone]
    }

print(clone_phylogeny_all_genotypes)


    print(clone_phylogeny_genotype)



}
