SIMULATOR_FULL_PHASE_2_copy_number_table <- function(package_sample){
    all_sample_genotype <- package_sample$all_sample_genotype
    all_sample_sampled_time <- package_sample$all_sample_sampled_time
    all_sample_ID <- package_sample$all_sample_ID
    sample_genotype_unique <- package_sample$sample_genotype_unique
    sample_genotype_unique_profile <- package_sample$sample_genotype_unique_profile
    #---Find the CN profiles for each cell in the sample
    sample_cell_ID <- c()
    sample_clone_ID <- all_sample_genotype
    sample_time <- all_sample_sampled_time

    for (i_cell in 1:length(all_sample_genotype)) {
        if (i_cell %% 10 == 0) {
            print(paste("Cell-", i_cell, sep = ""))
        }

        sample_ID <- all_sample_ID[i_cell]
        clone_ID <- sample_clone_ID[i_cell]
        i_clone <- which(sample_genotype_unique == clone_ID)[1]
        #       Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #       Add column for cell ID
        cell_ID <- paste(sample_ID, "-Library-", as.character(i_cell), "-", as.character(i_cell), sep = "")
        sample_cell_ID[i_cell] <- cell_ID
        cell_genotype_profile <- cbind(cell_genotype_profile, rep(cell_ID, nrow(cell_genotype_profile)))
        names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
        #       Update table of CN profiles for all cells in the sample
        if (i_cell == 1) {
            sample_genotype_profiles <- cell_genotype_profile
        } else {
            sample_genotype_profiles <- rbind(sample_genotype_profiles, cell_genotype_profile)
        }
    }

    package_sample$sample_genotype_profiles <- sample_genotype_profiles
    return(package_sample)
}
