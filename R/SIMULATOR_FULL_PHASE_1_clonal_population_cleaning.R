#' @export
SIMULATOR_FULL_PHASE_1_clonal_population_cleaning <- function(option_lite_memory) {
    #------------------------Function to find ancestors for given leaves
    find_ancestors <- function(clonal_ID_leaves) {
        clonal_ID_ancestors <- c()
        for (leaf in clonal_ID_leaves) {
            clone <- leaf
            while (clone > 0) {
                clonal_ID_ancestors <- c(clonal_ID_ancestors, clone)
                clone <- evolution_origin[clone]
            }
        }
        clonal_ID_ancestors <- unique(clonal_ID_ancestors)
        return(clonal_ID_ancestors)
    }
    #----------------------------------Find clonal populations to delete
    ind_clonal_population_delete <- which(clonal_population_next == 0)
    clonal_ID_delete <- clonal_ID_current[ind_clonal_population_delete]
    #--------------------------------Delete clonal populations & records
    if (length(ind_clonal_population_delete) != 0) {
        #---Delete clonal populations
        clonal_population_current <<- clonal_population_current[-ind_clonal_population_delete]
        clonal_population_next <<- clonal_population_next[-ind_clonal_population_delete]
        clonal_ID_current <<- clonal_ID_current[-ind_clonal_population_delete]
        #---Delete clonal records if requested
        if (option_lite_memory > 0) {
            if (option_lite_memory <= 2) {
                #---Only retain records for currently alive clones
                # genotype_list_ploidy_chrom
                genotype_list_ploidy_block[clonal_ID_delete] <<- -1
                genotype_list_ploidy_allele[clonal_ID_delete] <<- -1
                # genotype_list_WGD_count
                # genotype_list_driver_count
                # genotype_list_driver_map
                # genotype_list_selection_rate
                # genotype_list_DNA_length
                # genotype_list_prob_new_drivers
                # genotype_list_prob_CNAs
                # genotype_list_prob_CN_misseg_homolog
                # genotype_list_prob_CN_arm_misseg_homolog
            } else {
                #---Find ancestors of clones to delete
                clonal_ID_ancestors_current <- find_ancestors(clonal_ID_current)
                clonal_ID_ancestors_delete <- find_ancestors(clonal_ID_delete)
                clonal_ID_ancestors_delete <- clonal_ID_ancestors_delete[!(clonal_ID_ancestors_delete %in% clonal_ID_ancestors_current)]
                #---Only retain records for currently alive clones or their ancestors
                # genotype_list_ploidy_chrom
                genotype_list_ploidy_block[clonal_ID_ancestors_delete] <<- -1
                genotype_list_ploidy_allele[clonal_ID_ancestors_delete] <<- -1
                # genotype_list_WGD_count
                # genotype_list_driver_count
                # genotype_list_driver_map
                # genotype_list_selection_rate
                # genotype_list_DNA_length
                # genotype_list_prob_new_drivers
                # genotype_list_prob_CNAs
                # genotype_list_prob_CN_misseg_homolog
                # genotype_list_prob_CN_arm_misseg_homolog
            }
        }
    }
}
