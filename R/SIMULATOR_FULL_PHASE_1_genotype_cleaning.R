# ===================================DELETE NEW BUT UNNECESSARY GENOTYPES
#' @export
SIMULATOR_FULL_PHASE_1_genotype_cleaning <- function(genotype_to_react, genotype_daughter_1, genotype_daughter_2, position_to_react, position_daughter_1, position_daughter_2) {
    #-----------Find indicators of comparison between each pair of genotypes
    flag_mother_daughter1 <- SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_to_react, genotype_daughter_1)
    flag_mother_daughter2 <- SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_to_react, genotype_daughter_2)
    flag_daughter1_daughter2 <- SIMULATOR_FULL_PHASE_1_genotype_comparison(genotype_daughter_1, genotype_daughter_2)
    #--------------------------------------Delete new genotypes if necessary
    if ((flag_mother_daughter1 == 0) && (flag_mother_daughter2 == 0) && (flag_daughter1_daughter2 == 0)) {
        #       Case 1: all genotypes are different
        return
    } else {
        if ((flag_mother_daughter1 == 0) && (flag_mother_daughter2 == 1) && (flag_daughter1_daughter2 == 0)) {
            #       Case 2: combine genotype 2 into mother genotype
            N_clones <<- N_clones - 1

            evolution_origin <<- evolution_origin[-genotype_daughter_2]
            evolution_genotype_changes <<- evolution_genotype_changes[-genotype_daughter_2]

            genotype_list_ploidy_chrom <<- genotype_list_ploidy_chrom[-genotype_daughter_2]
            genotype_list_ploidy_allele <<- genotype_list_ploidy_allele[-genotype_daughter_2]
            genotype_list_ploidy_block <<- genotype_list_ploidy_block[-genotype_daughter_2]
            genotype_list_WGD_count <<- genotype_list_WGD_count[-genotype_daughter_2]
            genotype_list_driver_count <<- genotype_list_driver_count[-genotype_daughter_2]
            genotype_list_driver_map <<- genotype_list_driver_map[-genotype_daughter_2]
            genotype_list_DNA_length <<- genotype_list_DNA_length[-genotype_daughter_2]
            genotype_list_selection_rate <<- genotype_list_selection_rate[-genotype_daughter_2]
            genotype_list_prob_new_drivers <<- genotype_list_prob_new_drivers[-genotype_daughter_2]
            genotype_list_prob_CNAs <<- genotype_list_prob_CNAs[-genotype_daughter_2]
            if (mode_CN_misseg == "per_homolog") genotype_list_prob_CN_misseg_homolog <<- genotype_list_prob_CN_misseg_homolog[-genotype_daughter_2]
            if (mode_CN_arm_misseg == "per_homolog") genotype_list_prob_CN_arm_misseg_homolog <<- genotype_list_prob_CN_arm_misseg_homolog[-genotype_daughter_2]

            clonal_ID_current <<- clonal_ID_current[-position_daughter_2]
            clonal_population_current <<- clonal_population_current[-position_daughter_2]
            clonal_population_next <<- clonal_population_next[-position_daughter_2]

            if (genotype_daughter_1 > genotype_daughter_2) {
                genotype_daughter_1 <- genotype_daughter_1 - 1
            }
            genotype_daughter_2 <- genotype_to_react

            if (position_daughter_1 > position_daughter_2) {
                position_daughter_1 <- position_daughter_1 - 1
            }
            position_daughter_2 <- position_to_react
        } else {
            if ((flag_mother_daughter1 == 1) && (flag_mother_daughter2 == 0) && (flag_daughter1_daughter2 == 0)) {
                #       Case 3: combine genotype 1 into mother genotype
                N_clones <<- N_clones - 1

                evolution_origin <<- evolution_origin[-genotype_daughter_1]
                evolution_genotype_changes <<- evolution_genotype_changes[-genotype_daughter_1]

                genotype_list_ploidy_chrom <<- genotype_list_ploidy_chrom[-genotype_daughter_1]
                genotype_list_ploidy_allele <<- genotype_list_ploidy_allele[-genotype_daughter_1]
                genotype_list_ploidy_block <<- genotype_list_ploidy_block[-genotype_daughter_1]
                genotype_list_WGD_count <<- genotype_list_WGD_count[-genotype_daughter_1]
                genotype_list_driver_count <<- genotype_list_driver_count[-genotype_daughter_1]
                genotype_list_driver_map <<- genotype_list_driver_map[-genotype_daughter_1]
                genotype_list_DNA_length <<- genotype_list_DNA_length[-genotype_daughter_1]
                genotype_list_selection_rate <<- genotype_list_selection_rate[-genotype_daughter_1]
                genotype_list_prob_new_drivers <<- genotype_list_prob_new_drivers[-genotype_daughter_1]
                genotype_list_prob_CNAs <<- genotype_list_prob_CNAs[-genotype_daughter_1]
                if (mode_CN_misseg == "per_homolog") genotype_list_prob_CN_misseg_homolog <<- genotype_list_prob_CN_misseg_homolog[-genotype_daughter_1]
                if (mode_CN_arm_misseg == "per_homolog") genotype_list_prob_CN_arm_misseg_homolog <<- genotype_list_prob_CN_arm_misseg_homolog[-genotype_daughter_1]

                clonal_ID_current <<- clonal_ID_current[-position_daughter_1]
                clonal_population_current <<- clonal_population_current[-position_daughter_1]
                clonal_population_next <<- clonal_population_next[-position_daughter_1]

                if (genotype_daughter_2 > genotype_daughter_1) {
                    genotype_daughter_2 <- genotype_daughter_2 - 1
                }
                genotype_daughter_1 <- genotype_to_react

                if (position_daughter_2 > position_daughter_1) {
                    position_daughter_2 <- position_daughter_2 - 1
                    clonal_ID_current[position_daughter_2] <<- clonal_ID_current[position_daughter_2] - 1
                }
                position_daughter_1 <- position_to_react
            } else {
                if ((flag_mother_daughter1 == 0) && (flag_mother_daughter2 == 0) && (flag_daughter1_daughter2 == 1)) {
                    #       Case 4: combine genotype 2 into genotype 1
                    N_clones <<- N_clones - 1

                    evolution_origin <<- evolution_origin[-genotype_daughter_2]
                    evolution_genotype_changes <<- evolution_genotype_changes[-genotype_daughter_2]
                    genotype_list_ploidy_chrom <<- genotype_list_ploidy_chrom[-genotype_daughter_2]
                    genotype_list_ploidy_allele <<- genotype_list_ploidy_allele[-genotype_daughter_2]
                    genotype_list_ploidy_block <<- genotype_list_ploidy_block[-genotype_daughter_2]
                    genotype_list_WGD_count <<- genotype_list_WGD_count[-genotype_daughter_2]
                    genotype_list_driver_count <<- genotype_list_driver_count[-genotype_daughter_2]
                    genotype_list_driver_map <<- genotype_list_driver_map[-genotype_daughter_2]
                    genotype_list_DNA_length <<- genotype_list_DNA_length[-genotype_daughter_2]
                    genotype_list_selection_rate <<- genotype_list_selection_rate[-genotype_daughter_2]
                    genotype_list_prob_new_drivers <<- genotype_list_prob_new_drivers[-genotype_daughter_2]
                    genotype_list_prob_CNAs <<- genotype_list_prob_CNAs[-genotype_daughter_1]
                    if (mode_CN_misseg == "per_homolog") genotype_list_prob_CN_misseg_homolog <<- genotype_list_prob_CN_misseg_homolog[-genotype_daughter_1]
                    if (mode_CN_arm_misseg == "per_homolog") genotype_list_prob_CN_arm_misseg_homolog <<- genotype_list_prob_CN_arm_misseg_homolog[-genotype_daughter_1]

                    clonal_ID_current <<- clonal_ID_current[-position_daughter_2]
                    clonal_population_current <<- clonal_population_current[-position_daughter_2]
                    clonal_population_next <<- clonal_population_next[-position_daughter_2]

                    if (genotype_daughter_1 > genotype_daughter_2) {
                        genotype_daughter_1 <- genotype_daughter_1 - 1
                    }
                    genotype_daughter_2 <- genotype_daughter_1

                    if (position_daughter_1 > position_daughter_2) {
                        position_daughter_1 <- position_daughter_1 - 1
                    }
                    position_daughter_2 <- position_daughter_1
                } else {
                    if ((flag_mother_daughter1 == 1) && (flag_mother_daughter2 == 1) && (flag_daughter1_daughter2 == 1)) {
                        #       Case 5: combine both genotypes into mother genotype
                        N_clones <<- N_clones - 2

                        evolution_origin <<- evolution_origin[-c(genotype_daughter_1, genotype_daughter_2)]
                        evolution_genotype_changes <<- evolution_genotype_changes[-c(genotype_daughter_1, genotype_daughter_2)]

                        genotype_list_ploidy_chrom <<- genotype_list_ploidy_chrom[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_ploidy_allele <<- genotype_list_ploidy_allele[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_ploidy_block <<- genotype_list_ploidy_block[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_WGD_count <<- genotype_list_WGD_count[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_driver_count <<- genotype_list_driver_count[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_driver_map <<- genotype_list_driver_map[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_DNA_length <<- genotype_list_DNA_length[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_selection_rate <<- genotype_list_selection_rate[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_prob_new_drivers <<- genotype_list_prob_new_drivers[-c(genotype_daughter_1, genotype_daughter_2)]
                        genotype_list_prob_CNAs <<- genotype_list_prob_CNAs[-c(genotype_daughter_1, genotype_daughter_2)]
                        if (mode_CN_misseg == "per_homolog") genotype_list_prob_CN_misseg_homolog <<- genotype_list_prob_CN_misseg_homolog[-c(genotype_daughter_1, genotype_daughter_2)]
                        if (mode_CN_arm_misseg == "per_homolog") genotype_list_prob_CN_arm_misseg_homolog <<- genotype_list_prob_CN_arm_misseg_homolog[-c(genotype_daughter_1, genotype_daughter_2)]

                        clonal_ID_current <<- clonal_ID_current[-c(position_daughter_1, position_daughter_2)]
                        clonal_population_current <<- clonal_population_current[-c(position_daughter_1, position_daughter_2)]
                        clonal_population_next <<- clonal_population_next[-c(position_daughter_1, position_daughter_2)]

                        genotype_daughter_1 <- genotype_to_react
                        genotype_daughter_2 <- genotype_to_react

                        position_daughter_1 <- position_to_react
                        position_daughter_2 <- position_to_react
                    }
                }
            }
        }
    }
    output <- list()
    output[[1]] <- genotype_to_react
    output[[2]] <- genotype_daughter_1
    output[[3]] <- genotype_daughter_2
    output[[4]] <- position_to_react
    output[[5]] <- position_daughter_1
    output[[6]] <- position_daughter_2
    return(output)
}
