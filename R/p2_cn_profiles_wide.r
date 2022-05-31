p2_cn_profiles_wide <- function(simulation) {
    sample_cell_ID <- simulation$sample$sample_cell_ID

    all_sample_genotype <- simulation$sample$all_sample_genotype
    all_sample_sampled_time <- simulation$sample$all_sample_sampled_time
    all_sample_id <- simulation$sample$all_sample_ID
    sample_genotype_unique <- simulation$sample$sample_genotype_unique
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    #--------------Initiate the CN profiles with genome location indices
    cn_profiles_wide <- sample_genotype_unique_profile[[1]][, 1:3]
    #-------------------Find the CN profiles for each cell in the sample
    for (cell in 1:length(all_sample_genotype)) {
        clone_id <- all_sample_genotype[cell]
        i_clone <- which(sample_genotype_unique == clone_id)[1]
        #   Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #   Add column for cell ID
        cell_id <- sample_cell_ID[cell]
        cn_profiles_wide[cell_id] <- cell_genotype_profile$state
    }
    simulation$sample$cn_profiles_wide <- cn_profiles_wide
    return(simulation)
}
