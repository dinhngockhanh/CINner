#' @export
SIMULATOR_FULL_PHASE_1_clonal_population_cleaning <- function(option_lite_memory) {
    #----------------------------------Find clonal populations to delete
    vec_clonal_population_delete <- which(clonal_population_next == 0)
    #--------------------------------Delete clonal populations & records
    if (length(vec_clonal_population_delete) != 0) {
        #---Delete clonal records if requested
        if (option_lite_memory > 0) {
            vec_clonal_ID_delete <- clonal_ID_current[vec_clonal_population_delete]
            if (option_lite_memory <= 2) {
                #   Only retain records for currently alive clones
                # genotype_list_ploidy_chrom
                genotype_list_ploidy_block[vec_clonal_ID_delete] <<- -1
                genotype_list_ploidy_allele[vec_clonal_ID_delete] <<- -1
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
                stop("option_lite_memory must be 0, 1, or 2")
            }
        }
        #---Delete clonal populations
        clonal_population_current <<- clonal_population_current[-vec_clonal_population_delete]
        clonal_population_next <<- clonal_population_next[-vec_clonal_population_delete]
        clonal_ID_current <<- clonal_ID_current[-vec_clonal_population_delete]
    }
}
