#==============================PHASE 3: COPY-NUMBER PROFILES OF A SAMPLE
SIMULATOR_FULL_PHASE_3_main <- function(package_clonal_evolution,package_sample) {
#---------------------------------------------Input the clonal evolution
    T_current                                   <- package_clonal_evolution[[1]]
    N_clones                                    <- package_clonal_evolution[[4]]
    genotype_list_ploidy_chrom                  <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block                  <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele                 <- package_clonal_evolution[[7]]
    evolution_traj_time                         <- package_clonal_evolution[[13]]
    evolution_traj_divisions                    <- package_clonal_evolution[[14]]
    evolution_traj_clonal_ID                    <- package_clonal_evolution[[15]]
    evolution_traj_population                   <- package_clonal_evolution[[16]]
#-------------------------------------------------------Input the sample
    sample_cell_ID                              <- package_sample[[2]]
    sample_clone_ID                             <- package_sample[[3]]
#-----------------------------------Initialize phylogeny in hclust style
#   Initialize information to build phylogeny in hclust style
    hclust_row                                  <- 0
    hclust_nodes                                <- rep(0,1,2*N_sample-1)
    hclust_nodes[N_sample:2*N_sample-1]         <- (-1:-N_sample)
    hclust_labels                               <- sample_cell_ID
#   Initialize actual phylogeny in hclust style
    hclust_merge                                <- matrix(0,nrow=N_sample-1,ncol=2)
    hclust_height                               <- rep(0,1,N_sample-1)
#--------------------------------------Initialize phylogeny in our style
    phylogeny_origin                            <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_gens                      <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_genotypes                 <- vector("list",length=2*N_sample-1)
    phylogeny_genotype                          <- rep(0,length=2*N_sample-1)
    phylogeny_birthtime                         <- rep(0,length=2*N_sample-1)
    phylogeny_deathtime                         <- rep(0,length=2*N_sample-1)
#   Initialize the current list of node genotypes
    node_genotype_current                       <- sample_clone_ID
#   Initialize the current list of nodes in the sample phylogeny
    node_list_current                           <- N_sample:(2*N_sample-1)
#   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens[node_list_current]   <- 1
    for (node in N_sample:2*N_sample-1) {
        phylogeny_elapsed_genotypes[[node]]     <- node_genotype_current[node-N_sample+1]
    }
    phylogeny_genotype[node_list_current]       <- node_genotype_current
    phylogeny_deathtime[node_list_current]      <- T_current
#----------------------------------------Build the sample phylogeny tree
    for (i in seq(length(evolution_traj_divisions),1,-1)) {
#       Get time point
        time                        <- evolution_traj_time[i]

print(time)








    }







}
