#' @export
optimize_memory <- function(simulation) {
    # print(object.size(simulation))
    #--------------------------------Minimize clonal evolution variables
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    N_clones <- simulation$clonal_evolution$N_clones
    # genotype_list_ploidy_chrom<-simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    genotype_list_ploidy_allele <- simulation$clonal_evolution$genotype_list_ploidy_allele
    #   Find relevant clones to save
    if ("sample" %in% names(simulation)) {
        final_clones <- unique(simulation$sample$sample_clone_ID)
    } else {
        final_clones <- simulation$clonal_evolution$evolution_traj_clonal_ID[[length(simulation$clonal_evolution$evolution_traj_clonal_ID)]]
    }
    relevant_clones <- c()
    for (i in 1:length(final_clones)) {
        clone <- final_clones[i]
        while (clone > 0) {
            relevant_clones <- c(relevant_clones, clone)
            clone <- evolution_origin[clone]
        }
    }
    relevant_clones <- unique(relevant_clones)
    for (i in 1:N_clones) {
        if ((i %in% relevant_clones) == FALSE) {
            genotype_list_ploidy_block[[i]] <- -1
            genotype_list_ploidy_allele[[i]] <- -1
        }
    }
    # simulation$clonal_evolution$genotype_list_ploidy_chrom<-genotype_list_ploidy_chrom
    simulation$clonal_evolution$genotype_list_ploidy_block <- genotype_list_ploidy_block
    simulation$clonal_evolution$genotype_list_ploidy_allele <- genotype_list_ploidy_allele
    return(simulation)
}
