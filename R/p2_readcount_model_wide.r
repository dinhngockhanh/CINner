#' @export
p2_readcount_model_wide <- function(simulation, report_progress) {
    noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long

    sample_cell_ID <- simulation$sample$sample_cell_ID

    # all_sample_genotype <- simulation$sample$all_sample_genotype
    # all_sample_sampled_time <- simulation$sample$all_sample_sampled_time
    # all_sample_id <- simulation$sample$all_sample_ID
    # sample_genotype_unique <- simulation$sample$sample_genotype_unique
    # sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    #-------Initiate the readcount profiles with genome location indices
    noisy_cn_profiles_wide <- noisy_cn_profiles_long[which(noisy_cn_profiles_long$cell_id == sample_cell_ID[1]), 1:3]
    #------------Find the readcount profiles for each cell in the sample
    for (cell in 1:length(sample_cell_ID)) {
        cell_id <- sample_cell_ID[cell]
        #   Find the readcount profile for this cell
        cell_readcount_profile <- noisy_cn_profiles_long$reads[which(noisy_cn_profiles_long$cell_id == cell_id)]
        #   Add column for cell ID
        noisy_cn_profiles_wide[cell_id] <- cell_readcount_profile
    }
    #-------------------------------------------Output noisy CN profiles
    simulation$sample$noisy_cn_profiles_wide <- noisy_cn_profiles_wide
    return(simulation)
}
