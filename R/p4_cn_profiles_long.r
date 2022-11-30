#' @export
p4_cn_profiles_long <- function(simulation) {
    sample_cell_ID <- simulation$neutral_variations$sample$sample_cell_ID
    sample_clone_ID <- simulation$neutral_variations$sample$sample_clone_ID
    sample_genotype_unique <- simulation$neutral_variations$sample$sample_genotype_unique
    sample_genotype_unique_profile <- simulation$neutral_variations$sample$sample_genotype_unique_profile
    #-------------------Find the CN profiles for each cell in the sample
    cn_profiles_long_list <- vector("list", length = length(sample_clone_ID))
    for (cell in 1:length(sample_clone_ID)) {
        clone_ID <- sample_clone_ID[cell]
        i_clone <- which(sample_genotype_unique == clone_ID)[1]
        #   Find the CN profile for this cell
        cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
        #   Add column for cell ID
        cell_genotype_profile$cell_id <- sample_cell_ID[cell]
        names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
        #   Add record for this cell to list
        cn_profiles_long_list[[cell]] <- cell_genotype_profile
    }
    #-------------------------------Bind all cells' CN profiles together
    sample_genotype_profiles <- rbindlist(cn_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(sample_genotype_profiles) <- "data.frame"
    simulation$neutral_variations$sample$cn_profiles_long <- sample_genotype_profiles
    return(simulation)
}
