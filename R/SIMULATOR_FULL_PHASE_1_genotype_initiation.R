# ===========================INITIATE NEW GENOTYPES AFTER DRIVER/CN EVENT
#' @export
SIMULATOR_FULL_PHASE_1_genotype_initiation <- function(genotype_to_react) {
    #----------------------Create a new genotype for the first daughter cell
    N_clones <<- N_clones + 1
    genotype_daughter_1 <- N_clones

    position_daughter_1 <- length(clonal_population_current) + 1
    clonal_ID_current[position_daughter_1] <<- N_clones
    clonal_population_current[position_daughter_1] <<- 0
    clonal_population_next[position_daughter_1] <<- 0

    evolution_origin[N_clones] <<- genotype_to_react
    evolution_genotype_changes <<- c(evolution_genotype_changes, vector("list", length = 1))

    genotype_list_ploidy_chrom[[N_clones]] <<- genotype_list_ploidy_chrom[[genotype_to_react]]
    genotype_list_ploidy_allele[[N_clones]] <<- genotype_list_ploidy_allele[[genotype_to_react]]
    genotype_list_ploidy_block[[N_clones]] <<- genotype_list_ploidy_block[[genotype_to_react]]
    genotype_list_driver_count[N_clones] <<- genotype_list_driver_count[genotype_to_react]
    genotype_list_driver_map[[N_clones]] <<- genotype_list_driver_map[[genotype_to_react]]
    genotype_list_DNA_length[[N_clones]] <<- 0
    genotype_list_selection_rate[N_clones] <<- 0
    genotype_list_prob_new_drivers[N_clones] <<- 0
    #---------------------Create a new genotype for the second daughter cell
    N_clones <<- N_clones + 1
    genotype_daughter_2 <- N_clones

    position_daughter_2 <- length(clonal_population_current) + 1
    clonal_ID_current[position_daughter_2] <<- N_clones
    clonal_population_current[position_daughter_2] <<- 0
    clonal_population_next[position_daughter_2] <<- 0

    evolution_origin[N_clones] <<- genotype_to_react
    evolution_genotype_changes <<- c(evolution_genotype_changes, vector("list", length = 1))

    genotype_list_ploidy_chrom[[N_clones]] <<- genotype_list_ploidy_chrom[[genotype_to_react]]
    genotype_list_ploidy_allele[[N_clones]] <<- genotype_list_ploidy_allele[[genotype_to_react]]
    genotype_list_ploidy_block[[N_clones]] <<- genotype_list_ploidy_block[[genotype_to_react]]
    genotype_list_driver_count[N_clones] <<- genotype_list_driver_count[genotype_to_react]
    genotype_list_driver_map[[N_clones]] <<- genotype_list_driver_map[[genotype_to_react]]
    genotype_list_DNA_length[[N_clones]] <<- 0
    genotype_list_selection_rate[N_clones] <<- 0
    genotype_list_prob_new_drivers[N_clones] <<- 0
    #-------------------------------------Output the new initiated genotypes
    output <<- list()
    output[[1]] <<- genotype_daughter_1
    output[[2]] <<- genotype_daughter_2
    output[[3]] <<- position_daughter_1
    output[[4]] <<- position_daughter_2
    return(output)
}
