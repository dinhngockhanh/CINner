#' @export
SIMULATOR_FULL_PHASE_4_main <- function(package_clonal_evolution, package_sample, package_sample_phylogeny, report_progress) {
    #-----------------------------------------------Input CN event rates
    prob_CN_WGD <- prob_neutral_CN_whole_genome_duplication
    prob_CN_misseg <- prob_neutral_CN_missegregation
    prob_CN_arm_misseg <- prob_neutral_CN_chrom_arm_missegregation
    prob_CN_foc_amp <- prob_neutral_CN_focal_amplification
    prob_CN_foc_del <- prob_neutral_CN_focal_deletion
    prob_CN_cnloh_i <- prob_neutral_CN_cnloh_interstitial
    prob_CN_cnloh_t <- prob_neutral_CN_cnloh_terminal
    #---------------------------------------------------Input the sample
    sample_cell_ID <- package_sample$sample_cell_ID
    # all_sample_genotype <- package_sample$all_sample_genotype
    sample_clone_ID <- package_sample$sample_clone_ID
    sample_clone_ID_letters <- package_sample$sample_clone_ID_letters
    table_clone_ID_vs_letters <- package_sample$table_clone_ID_vs_letters
    sample_time <- package_sample$sample_time
    #-----------------------------------------Input the sample phylogeny
    phylogeny_origin <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_origin
    phylogeny_elapsed_gens <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_elapsed_gens
    phylogeny_elapsed_genotypes <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_elapsed_genotypes
    phylogeny_genotype <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_genotype
    phylogeny_birthtime <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_birthtime
    phylogeny_deathtime <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_deathtime
    phylogeny_order <- package_sample_phylogeny$package_cell_phylogeny$phylogeny_order
    #-----------------------Initialize phylogeny with neutral variations
    neuvar_phylogeny_origin <- phylogeny_origin
    neuvar_phylogeny_elapsed_gens <- phylogeny_elapsed_gens
    neuvar_phylogeny_elapsed_genotypes <- vector("list", length = length(phylogeny_elapsed_genotypes))
    neuvar_phylogeny_genotype <- rep(0, length = length(phylogeny_genotype))
    neuvar_phylogeny_birthtime <- phylogeny_birthtime
    neuvar_phylogeny_deathtime <- phylogeny_deathtime
    neuvar_phylogeny_order <- phylogeny_order
    #--------------Get list of future genotypes growing from each branch
    list_progeny_genotypes <- vector("list", length = length(phylogeny_genotype))
    for (node in seq(length(phylogeny_genotype), 1, -1)) {
        list_progeny_genotypes[[node]] <- unique(c(list_progeny_genotypes[[node]], phylogeny_elapsed_genotypes[[node]]))
        node_mother <- phylogeny_origin[node]
        if (node_mother > 0) {
            list_progeny_genotypes[[node_mother]] <- unique(c(list_progeny_genotypes[[node_mother]], list_progeny_genotypes[[node]]))
        }
    }
    #--------------------------Simulate the subclonal neutral variations
    if (report_progress == TRUE) {
        pb <- txtProgressBar(
            min = 1, max = length(neuvar_phylogeny_origin),
            style = 3, width = 50, char = "="
        )
    }
    for (branch in 1:length(neuvar_phylogeny_origin)) {
        if (report_progress == TRUE) {
            setTxtProgressBar(pb, branch)
        }
        # ==Get mother genotype with neutral variations
        neuvar_mother_branch <- neuvar_phylogeny_origin[branch]
        if (neuvar_mother_branch <= 0) {
            neuvar_mother_genotype <- phylogeny_elapsed_genotypes[[branch]][1]
        } else {
            neuvar_mother_genotype <- neuvar_phylogeny_genotype[neuvar_mother_branch]
        }
        # ==Get list of all future genotypes growing from this branch
        progeny_genotypes <- list_progeny_genotypes[[branch]]
        # ==Get list of elapsed genotypes in original phylogeny
        original_elapsed_genotypes <- phylogeny_elapsed_genotypes[[branch]]
        # ==Create new elapsed genotypes with neutral variations
        neuvar_elapsed_genotypes <- c()
        genotype_current <- neuvar_mother_genotype
        for (elapsed_gen in 1:length(original_elapsed_genotypes)) {
            #---Get list of chromosomes to exclude from neutral variations
            if (elapsed_gen == length(original_elapsed_genotypes)) {
                future_genotypes <- progeny_genotypes
            } else {
                future_genotypes <- setdiff(original_elapsed_genotypes[(elapsed_gen + 1):length(original_elapsed_genotypes)], original_elapsed_genotypes[elapsed_gen])
                future_genotypes <- unique(c(future_genotypes, progeny_genotypes))
            }
            chroms_excluded <- c()
            if (length(future_genotypes) > 0) {
                for (i in 1:length(future_genotypes)) {
                    genotype <- future_genotypes[i]
                    events <- evolution_genotype_changes[[genotype]]
                    if (length(events) == 0) {
                        next
                    }
                    for (j in 1:length(events)) {
                        # if ((events[[j]][1] == "missegregation") & (events[[j]][4] == -1)) {
                        #     chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        # } else if ((events[[j]][1] == "chromosome-arm-missegregation") & (events[[j]][5] == -1)) {
                        #     chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        # } else if (events[[j]][1] == "focal-deletion") {
                        #     chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        # } else if (events[[j]][1] == "cnloh-terminal") {
                        #     chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        # } else if (events[[j]][1] == "cnloh-interstitial") {
                        #     chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        # }



                        if (events[[j]][1] != "whole-genome-duplication") {
                            chroms_excluded <- c(chroms_excluded, strtoi(events[[j]][2]))
                        }
                    }
                }
            }
            #---Perform the pre-determined events (if any)
            if (elapsed_gen > 1) {
                previous_genotype <- original_elapsed_genotypes[elapsed_gen - 1]
            } else if (branch > 1) {
                elapsed_gens_mother <- phylogeny_elapsed_genotypes[[neuvar_phylogeny_origin[branch]]]
                previous_genotype <- elapsed_gens_mother[length(elapsed_gens_mother)]
            } else {
                previous_genotype <- original_elapsed_genotypes[elapsed_gen]
            }
            if (original_elapsed_genotypes[elapsed_gen] != previous_genotype) {
                #   Create a new genotype
                genotype_parent <- genotype_current
                output <- SIMULATOR_FULL_PHASE_1_genotype_initiation(genotype_current, n_daughters = 1)
                genotype_current <- output[[1]]
                #   Perform driver and CN events
                events <- evolution_genotype_changes[[original_elapsed_genotypes[elapsed_gen]]]
                for (j in 1:length(events)) {
                    event <- events[[j]]
                    if (event[1] == "new-driver") {
                        SIMULATOR_FULL_PHASE_1_drivers(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "whole-genome-duplication") {
                        SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current
                        )
                    } else if (event[1] == "missegregation") {
                        SIMULATOR_FULL_PHASE_1_CN_missegregation(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "chromosome-arm-missegregation") {
                        SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "focal-amplification") {
                        SIMULATOR_FULL_PHASE_1_CN_focal_amplification(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "focal-deletion") {
                        SIMULATOR_FULL_PHASE_1_CN_focal_deletion(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "cnloh-interstitial") {
                        SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            event = event
                        )
                    } else if (event[1] == "cnloh-terminal") {
                        SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            event = event
                        )
                    }
                }
            }
            selection_rate <- SIMULATOR_FULL_PHASE_1_selection_rate(
                genotype_list_driver_count[genotype_current],
                genotype_list_driver_map[[genotype_current]],
                genotype_list_ploidy_chrom[[genotype_current]],
                genotype_list_ploidy_block[[genotype_current]],
                genotype_list_ploidy_allele[[genotype_current]]
            )
            genotype_list_selection_rate[genotype_current] <<- selection_rate
            #---Perform neutral variations
            genotype_current_tmp <- genotype_current
            selection_rate_tmp <- genotype_list_selection_rate[genotype_current]
            selection_rate <- 0
            while ((selection_rate <= 0) & (selection_rate < selection_rate_tmp)) {
                genotype_current <- genotype_current_tmp
                flag_whole_genome_duplication <- as.numeric(runif(1) < prob_CN_WGD)
                flag_missegregation <- as.numeric(runif(1) < prob_CN_misseg)
                flag_chrom_arm_missegregation <- as.numeric(runif(1) < prob_CN_arm_misseg)
                flag_amplification <- as.numeric(runif(1) < prob_CN_foc_amp)
                flag_deletion <- as.numeric(runif(1) < prob_CN_foc_del)
                flag_cnloh_interstitial <- as.numeric(runif(1) < prob_CN_cnloh_i)
                flag_cnloh_terminal <- as.numeric(runif(1) < prob_CN_cnloh_t)
                vec_flag <- c(flag_whole_genome_duplication, flag_missegregation, flag_chrom_arm_missegregation, flag_amplification, flag_deletion, flag_cnloh_interstitial, flag_cnloh_terminal)
                if (max(vec_flag) == 1) {
                    #   Create a new genotype
                    genotype_parent <- genotype_current
                    output <- SIMULATOR_FULL_PHASE_1_genotype_initiation(genotype_current, n_daughters = 1)
                    genotype_current <- output[[1]]
                    #   Perform neutral CN events
                    if (flag_whole_genome_duplication == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current
                        )
                    }
                    if (flag_missegregation == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_missegregation(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                    if (flag_chrom_arm_missegregation == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation(
                            genotype_to_react = genotype_parent,
                            genotype_daughter_1 = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                    if (flag_amplification == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_focal_amplification(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                    if (flag_deletion == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_focal_deletion(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                    if (flag_cnloh_interstitial == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                    if (flag_cnloh_terminal == 1) {
                        SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal(
                            genotype_to_react = genotype_parent,
                            genotype_daughter = genotype_current,
                            chromosomes_excluded = chroms_excluded
                        )
                    }
                }
                #   Compute the selection rate for the new genotype
                selection_rate <- SIMULATOR_FULL_PHASE_1_selection_rate(
                    genotype_list_driver_count[genotype_current],
                    genotype_list_driver_map[[genotype_current]],
                    genotype_list_ploidy_chrom[[genotype_current]],
                    genotype_list_ploidy_block[[genotype_current]],
                    genotype_list_ploidy_allele[[genotype_current]]
                )
                genotype_list_selection_rate[genotype_current] <<- selection_rate
            }
            #---Update elapsed genotypes
            neuvar_elapsed_genotypes <- c(neuvar_elapsed_genotypes, genotype_current)
        }
        # ==Update phylogeny records with neutral variations
        neuvar_phylogeny_elapsed_genotypes[[branch]] <- neuvar_elapsed_genotypes
        neuvar_phylogeny_genotype[branch] <- genotype_current
    }
    if (report_progress == TRUE) {
        cat("\n")
    }
    #------------Update clonal evolution package with new neutral clones
    package_clonal_evolution$N_clones <- N_clones
    package_clonal_evolution$genotype_list_ploidy_chrom <- genotype_list_ploidy_chrom
    package_clonal_evolution$genotype_list_ploidy_allele <- genotype_list_ploidy_allele
    package_clonal_evolution$genotype_list_ploidy_block <- genotype_list_ploidy_block
    package_clonal_evolution$genotype_list_driver_count <- genotype_list_driver_count
    package_clonal_evolution$genotype_list_driver_map <- genotype_list_driver_map
    package_clonal_evolution$genotype_list_selection_rate <- genotype_list_selection_rate
    package_clonal_evolution$evolution_origin <- evolution_origin
    package_clonal_evolution$evolution_genotype_changes <- evolution_genotype_changes
    #------------------------------Get the CN profiles for sampled cells
    neuvar_sample_cell_ID <- sample_cell_ID
    neuvar_sample_clone_ID <- neuvar_phylogeny_genotype[(length(neuvar_phylogeny_genotype) - length(sample_cell_ID) + 1):length(neuvar_phylogeny_genotype)]
    neuvar_sample_genotype_unique <- unique(neuvar_sample_clone_ID)
    neuvar_sample_genotype_unique_profile <- list()
    for (i_clone in 1:length(neuvar_sample_genotype_unique)) {
        clone_ID <- neuvar_sample_genotype_unique[i_clone]
        genotype_unique_profile <- get_cn_profile(package_clonal_evolution, clone_ID)
        neuvar_sample_genotype_unique_profile[[i_clone]] <- genotype_unique_profile
    }
    #-----------------------------Output package of data from simulation
    cell_phylogeny <- list()
    cell_phylogeny$neuvar_phylogeny_origin <- neuvar_phylogeny_origin
    cell_phylogeny$neuvar_phylogeny_elapsed_gens <- neuvar_phylogeny_elapsed_gens
    cell_phylogeny$neuvar_phylogeny_elapsed_genotypes <- neuvar_phylogeny_elapsed_genotypes
    cell_phylogeny$neuvar_phylogeny_genotype <- neuvar_phylogeny_genotype
    cell_phylogeny$neuvar_phylogeny_birthtime <- neuvar_phylogeny_birthtime
    cell_phylogeny$neuvar_phylogeny_deathtime <- neuvar_phylogeny_deathtime
    cell_phylogeny$neuvar_phylogeny_order <- neuvar_phylogeny_order

    sample <- list()
    sample$sample_cell_ID <- neuvar_sample_cell_ID
    sample$sample_clone_ID <- neuvar_sample_clone_ID
    sample$sample_genotype_unique <- neuvar_sample_genotype_unique
    sample$sample_genotype_unique_profile <- neuvar_sample_genotype_unique_profile

    output <- list()
    output$cell_phylogeny <- cell_phylogeny
    output$sample <- sample

    return(output)
}
