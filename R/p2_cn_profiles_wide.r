p2_cn_profiles_wide <- function(simulation) {
    all_sample_genotype <- simulation$sample$all_sample_genotype
    all_sample_sampled_time <- simulation$sample$all_sample_sampled_time
    all_sample_id <- simulation$sample$all_sample_ID
    sample_genotype_unique <- simulation$sample$sample_genotype_unique
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile

    print(sample_genotype_unique_profile[[1]])
    #--------------Initiate the CN profiles with genome location indices
    cn_profiles_wide <- sample_genotype_unique_profile[[1]][, 1:3]
    #-------------------Find the CN profiles for each cell in the sample
    pb <- txtProgressBar(
        min = 0, max = length(all_sample_genotype),
        style = 3, width = 50, char = "="
    )
    for (cell in 1:length(all_sample_genotype)) {
        setTxtProgressBar(pb, cell)
        sample_id <- all_sample_id[cell]
        clone_id <- all_sample_genotype[cell]
        i_clone <- which(sample_genotype_unique == clone_id)[1]
        #       Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #       Add column for cell ID
        cell_id <- paste(sample_id, "-Library-", as.character(cell), "-", as.character(cell), sep = "")

        cn_profiles_wide[cell_id] <- cell_genotype_profile$state
    }
    cat("\n")
    simulation$sample$cn_profiles_wide <- cn_profiles_wide
    return(simulation)
}
