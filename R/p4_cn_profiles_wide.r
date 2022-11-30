#' @export
p4_cn_profiles_wide <- function(simulation) {
    sample_cell_ID <- simulation$neutral_variations$sample$sample_cell_ID
    sample_clone_ID <- simulation$neutral_variations$sample$sample_clone_ID
    sample_genotype_unique <- simulation$neutral_variations$sample$sample_genotype_unique
    sample_genotype_unique_profile <- simulation$neutral_variations$sample$sample_genotype_unique_profile
    #--------------Initiate the CN profiles with genome location indices
    cn_profiles_wide <- sample_genotype_unique_profile[[1]][, 1:3]
    #-------------------Find the CN profiles for each cell in the sample
    for (cell in 1:length(sample_clone_ID)) {
        clone_id <- sample_clone_ID[cell]
        i_clone <- which(sample_genotype_unique == clone_id)[1]
        #   Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #   Add column for cell ID
        cell_id <- sample_cell_ID[cell]
        cn_profiles_wide[cell_id] <- cell_genotype_profile$state
    }
    simulation$neutral_variations$sample$cn_profiles_wide <- cn_profiles_wide
    return(simulation)
}
