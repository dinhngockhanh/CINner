p2_cn_profiles_long <- function(simulation) {
    # all_sample_genotype <- simulation$sample$all_sample_genotype
    # all_sample_sampled_time <- simulation$sample$all_sample_sampled_time
    # all_sample_ID <- simulation$sample$all_sample_ID
    # sample_genotype_unique <- simulation$sample$sample_genotype_unique
    # sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    # #-------------------Find the CN profiles for each cell in the sample
    # sample_cell_ID <- c()
    # sample_clone_ID <- all_sample_genotype
    # sample_time <- all_sample_sampled_time
    # pb <- txtProgressBar(
    #     min = 0, max = length(all_sample_genotype),
    #     style = 3, width = 50, char = "="
    # )
    # #   Make list of CN profiles for each cell
    # cn_profiles_long_list<-vector("list", length = length(all_sample_genotype))
    # for (cell in 1:length(all_sample_genotype)) {
    #     setTxtProgressBar(pb, cell)
    #     sample_ID <- all_sample_ID[cell]
    #     clone_ID <- sample_clone_ID[cell]
    #     i_clone <- which(sample_genotype_unique == clone_ID)[1]
    #     #       Find the CN profile for this cell
    #     cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
    #     #       Add column for cell ID
    #     cell_ID <- paste(sample_ID, "-Library-", as.character(cell), "-", as.character(cell), sep = "")
    #     sample_cell_ID[cell] <- cell_ID
    #     cell_genotype_profile <- cbind(cell_genotype_profile, rep(cell_ID, nrow(cell_genotype_profile)))
    #     names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
    #     #       Update table of CN profiles for all cells in the sample
    #
    #     cn_profiles_long_list[[cell]] <- cell_genotype_profile
    #
    #     if (cell == 1) {
    #         cn_profiles_long <- cell_genotype_profile
    #     } else {
    #         cn_profiles_long <- rbind(cn_profiles_long, cell_genotype_profile)
    #     }
    # }
    # cat("\n")

    cn_profiles_long <- simulation$sample$sample_genotype_profiles
    class(cn_profiles_long) <- "data.frame"
    simulation$sample$cn_profiles_long <- cn_profiles_long
    return(simulation)
}
