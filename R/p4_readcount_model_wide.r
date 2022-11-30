#' @export
p4_readcount_model_wide <- function(simulation, report_progress) {
    noisy_cn_profiles_long <- simulation$neutral_variations$sample$noisy_cn_profiles_long
    sample_cell_ID <- simulation$neutral_variations$sample$sample_cell_ID
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
    simulation$neutral_variations$sample$noisy_cn_profiles_wide <- noisy_cn_profiles_wide
    return(simulation)
}
