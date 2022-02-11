# #----------------------------------------Initialize the phylogeny record
#     phylogeny_origin                <- rep(0,length=2*N_sample-1)
#     phylogeny_elapsed_gens          <- rep(0,length=2*N_sample-1)
#     phylogeny_elapsed_genotypes     <- vector("list",length=2*N_sample-1)
#     phylogeny_genotype              <- rep(0,length=2*N_sample-1)
#     phylogeny_birthtime             <- rep(0,length=2*N_sample-1)
#     phylogeny_deathtime             <- rep(0,length=2*N_sample-1)
# #-------------------------------Find a random sample of final population
# #   Initialize the current list of nodes in the sample phylogeny
#     node_list_current               <- N_sample:(2*N_sample-1)
# #   Initialize data for leaves of sample phylogeny
#     phylogeny_elapsed_gens[node_list_current]   <- 1
#     for (node in N_sample:2*N_sample-1) {
#         phylogeny_elapsed_genotypes[[node]]     <- node_genotype_current[node-N_sample+1]
#     }
#     phylogeny_genotype[node_list_current]       <- node_genotype_current
#     phylogeny_deathtime[node_list_current]      <- T_current
