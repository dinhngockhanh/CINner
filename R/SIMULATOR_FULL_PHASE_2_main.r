#==============================================PHASE 2: SAMPLE PHYLOGENY
SIMULATOR_FULL_PHASE_2_main <- function(package_clonal_evolution) {
#---------------------------------------------Input the clonal evolution
    T_current                       <- package_clonal_evolution[[1]]
    N_clones                        <- package_clonal_evolution[[4]]
    genotype_list_ploidy_chrom      <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block      <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele     <- package_clonal_evolution[[7]]
    evolution_traj_time             <- package_clonal_evolution[[13]]
    evolution_traj_divisions        <- package_clonal_evolution[[14]]
    evolution_traj_clonal_ID        <- package_clonal_evolution[[15]]
    evolution_traj_population       <- package_clonal_evolution[[16]]
#----------------------------------------Initialize the phylogeny record
    phylogeny_origin                <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_gens          <- rep(0,length=2*N_sample-1)
    phylogeny_elapsed_genotypes     <- vector("list",length=2*N_sample-1)
    phylogeny_genotype              <- rep(0,length=2*N_sample-1)
    phylogeny_birthtime             <- rep(0,length=2*N_sample-1)
    phylogeny_deathtime             <- rep(0,length=2*N_sample-1)
#-------------------------------Find a random sample of final population
#   Initialize the current list of nodes in the sample phylogeny
    node_list_current               <- N_sample:(2*N_sample-1)
#   Initialize the current list of genoytpes of the nodes
    final_clonal_ID                 <- tail(evolution_traj_clonal_ID,1)[[1]]
    final_clonal_population         <- tail(evolution_traj_population,1)[[1]]
    final_population                <- c()
    for (i in 1:length(final_clonal_ID)) {
        clone                       <- final_clonal_ID[i]
        clonal_population           <- final_clonal_population[i]
        final_population            <- c(final_population,rep(clone,1,clonal_population)) # append clone to final_population
    }
    node_genotype_current           <- sample(x=final_population,size=N_sample,replace=FALSE)
#   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens[node_list_current]   <- 1
    for (node in N_sample:2*N_sample-1) {
        phylogeny_elapsed_genotypes[[node]]     <- node_genotype_current[node-N_sample+1]
    }
    phylogeny_genotype[node_list_current]       <- node_genotype_current
    phylogeny_deathtime[node_list_current]      <- T_current





#---------------------------------Create CN object for the sampled cells
    sample_genotype                 <- phylogeny_genotype[N_sample:(2*N_sample-1)]
#---Find the CN profiles for each clone found in the sample
    sample_genotype_unique          <- unique(sample_genotype)
    sample_genotype_unique_profile  <- list()
    for(i_clone in 1:1){
    # for(i_clone in 1:length(sample_genotype_unique)){
#       Extract CN information for the clone from clonal evolution data
        clone_ID                    <- sample_genotype_unique[i_clone]
        ploidy_chrom                <- genotype_list_ploidy_chrom[[clone_ID]]
        ploidy_block                <- genotype_list_ploidy_block[[clone_ID]]
        ploidy_allele               <- genotype_list_ploidy_allele[[clone_ID]]
#       Build the CN profile in SIGNALS style for the clone
        vec_clone_chr               <- c()
        vec_clone_start             <- c()
        vec_clone_end               <- c()
        vec_clone_copy              <- c()
        vec_clone_state             <- c()
        vec_clone_Min               <- c()
        vec_clone_Maj               <- c()
        for (chrom in 1:1){
        # for (chrom in 1:N_chromosomes){
            chrom_block_count       <- vec_CN_block_no[chrom]
            chrom_ploidy            <- ploidy_chrom[chrom]
#           Find location information of each chromosome block
            vec_chr                 <- rep(chrom,1,chrom_block_count)
            vec_start               <- seq(0,size_CN_block_DNA*(chrom_block_count-1),by=size_CN_block_DNA)+1
            vec_end                 <- seq(size_CN_block_DNA,size_CN_block_DNA*chrom_block_count,by=size_CN_block_DNA)





print(vec_chr)
print(vec_start)
print(vec_end)

        }


# print('---------------------------------------------------------------')
# print(clone_ID)
# print(ploidy_chrom)
# print(ploidy_block[[1]])
# print(ploidy_allele[[1]])

    }




# print(sample_genotype)
# print(sample_genotype_unique)


# print(N_sample)
# print(genotype_list_ploidy_block[[1]])
# print(genotype_list_ploidy_allele[[1]])


}
