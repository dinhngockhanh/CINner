SIMULATOR_FULL_PHASE_2_copy_number_table <- function(package_sample) {
    all_sample_genotype <- package_sample$all_sample_genotype
    all_sample_sampled_time <- package_sample$all_sample_sampled_time
    all_sample_ID <- package_sample$all_sample_ID
    sample_genotype_unique <- package_sample$sample_genotype_unique
    sample_genotype_unique_profile <- package_sample$sample_genotype_unique_profile
    #---Find the CN profiles for each cell in the sample
    sample_cell_ID <- c()
    sample_clone_ID <- all_sample_genotype
    sample_time <- all_sample_sampled_time

    pb <- txtProgressBar(
        min = 0,
        max = length(all_sample_genotype),
        style = 3,
        width = 50,
        char = "="
    )
    #   Make list of CN profiles for each cell
    sample_genotype_profiles_list <- vector("list", length = length(all_sample_genotype))
    for (i_cell in 1:length(all_sample_genotype)) {
        sample_ID <- all_sample_ID[i_cell]
        clone_ID <- sample_clone_ID[i_cell]
        i_clone <- which(sample_genotype_unique == clone_ID)[1]
        #       Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #       Add column for cell ID
        cell_ID <- paste(sample_ID, "-Library-", as.character(i_cell), "-", as.character(i_cell), sep = "")
        sample_cell_ID[i_cell] <- cell_ID
        cell_genotype_profile$cell_id <- cell_ID
        names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
        #       Add record for this cell to list
        sample_genotype_profiles_list[[i_cell]] <- cell_genotype_profile
        setTxtProgressBar(pb, i_cell)
    }
    cat("\n")

    #   Bind all cells' CN profiles together
    sample_genotype_profiles <- rbindlist(sample_genotype_profiles_list, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(sample_genotype_profiles) <- "data.frame"
    package_sample$sample_genotype_profiles <- sample_genotype_profiles
    return(package_sample)
}
