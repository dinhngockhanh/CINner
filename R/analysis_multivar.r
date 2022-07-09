#' @export
simulator_multivar <- function(model_prefix = "",
                               model_variables_base = list(),
                               var1_name = "",
                               var1_vals = c(),
                               var2_name = "",
                               var2_vals = c(),
                               n_simulations_per_batch = 0,
                               compute_parallel = FALSE,
                               stage_final = 3,
                               n_clones_min = 1,
                               n_clones_max = Inf,
                               seed = Inf) {
    # ================================================CREATE SIMULATIONS
    pb <- txtProgressBar(
        min = 1, max = length(var1_vals) * length(var2_vals),
        style = 3, width = 50, char = "="
    )
    # for (row in 1:length(var1_vals)) {
    #     for (col in 1:length(var2_vals)) {
    for (row in length(var1_vals):1) {
        for (col in length(var2_vals):1) {
            setTxtProgressBar(pb, ((row - 1) * length(var2_vals) + col))
            #   Initialize parameter set
            model_variables <- model_variables_base
            model_name <- paste(model_prefix, "-", row, "-", col, sep = "")
            #   Fix variable 1 in parameter set
            var1 <- var1_vals[row]
            if (var1_name == "rate_WGD") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var1
            } else {
                if (var1_name == "rate_missegregation") {
                    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- var1
                }
            }
            #   Fix variable 2 in parameter set
            var2 <- var2_vals[col]
            if (var2_name == "rate_WGD") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var2
            } else {
                if (var2_name == "rate_missegregation") {
                    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- var2
                }
            }
            #   Save parameter set
            model_name <- paste(model_prefix, "-", row, "-", col, sep = "")
            SAVE_model_variables(model_name = model_name, model_variables = model_variables)
            #   Run simulations for the entry
            n_simulations <- n_simulations_per_batch
            if (seed == Inf) {
                set.seed(Sys.time())
            } else {
                set.seed(seed)
            }
            simulator_full_program(
                model = model_name,
                n_simulations = n_simulations,
                stage_final = stage_final,
                n_clones_min = n_clones_min,
                n_clones_max = n_clones_max,
                save_simulation = TRUE,
                # save_simulation = FALSE,
                save_newick_tree = FALSE,
                save_cn_profile = FALSE,
                # format_cn_profile = "both",
                internal_nodes_cn_info = FALSE,
                model_readcount = FALSE,
                apply_HMM = FALSE,
                report_progress = FALSE,
                compute_parallel = compute_parallel,
                pseudo_corrected_readcount = FALSE,
                seed = seed
            )
        }
    }
    cat("\n")
}
#' @export
statistics_multivar <- function(model_prefix = "",
                                var1_name = "",
                                var1_vals = c(),
                                var2_name = "",
                                var2_vals = c(),
                                n_simulations_per_batch = 0,
                                stage_final = 3) {
    # =====================================SIMULATE MATRIX OF STATISTICS
    mat_simulation_statistics <- list()
    pb <- txtProgressBar(
        min = 1, max = length(var1_vals) * length(var2_vals),
        style = 3, width = 50, char = "="
    )
    for (row in 1:length(var1_vals)) {
        mat_simulation_statistics[[row]] <- list()
        for (col in 1:length(var2_vals)) {
            setTxtProgressBar(pb, ((row - 1) * length(var2_vals) + col))
            model <- paste(model_prefix, "-", row, "-", col, sep = "")
            SIMULATOR_VARIABLES_for_simulation(model)
            #   Find statistics from each simulation
            batch_statistics <- list()
            for (i in 1:n_simulations_per_batch) {
                filename <- paste(model, "_simulation_", i, ".rda", sep = "")
                load(filename)
                simulation <- simulation_stats(simulation, stage_final)
                batch_statistics[[i]] <- simulation$statistics
            }
            mat_simulation_statistics[[row]][[col]] <- batch_statistics
        }
    }
    cat("\n")
    # =========================================SAVE MATRIX OF STATISTICS
    save(mat_simulation_statistics, file = paste(model_prefix, "_statistics.rda", sep = ""))
}
#' @export
simulation_stats <- function(simulation, stage_final) {
    simulation$statistics <- list()
    #-------------------------Statistics from phase 1 (clonal evolution)
    if (stage_final >= 1) {
        #---True clone counts in whole population
        true_n_clones_whole_pop_times <- simulation$clonal_evolution$evolution_traj_time
        true_n_clones_whole_pop_values <- rep(0, length = length(true_n_clones_whole_pop_times))
        for (step in 1:length(true_n_clones_whole_pop_values)) {
            true_n_clones_whole_pop_values[step] <- length(simulation$clonal_evolution$evolution_traj_clonal_ID[[step]])
        }
        simulation$statistics$true_n_clones_whole_pop_times <- true_n_clones_whole_pop_times
        simulation$statistics$true_n_clones_whole_pop_values <- true_n_clones_whole_pop_values
    }
    #---------------------------------Statistics from phase 2 (sampling)
    if (stage_final >= 2) {
        #---Prepare dataframe for statistics
        stats_sample <- Table_sampling
        #---True clone counts in sample
        stats_sample$true_n_clones <- 0
        for (sample in 1:nrow(stats_sample)) {
            vec_loc <- which(simulation$sample$all_sample_ID == stats_sample$Sample_ID[sample])
            stats_sample$true_n_clones[sample] <- length(unique(simulation$sample$all_sample_genotype[vec_loc]))
        }
        #---Total-CN-based clone counts in sample
        stats_sample$total_cn_based_n_clones <- 0
        genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
        genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
        for (sample in 1:nrow(stats_sample)) {
            vec_loc <- which(simulation$sample$all_sample_ID == stats_sample$Sample_ID[sample])
            vec_unique_genotypes <- unique(simulation$sample$all_sample_genotype[vec_loc])
            #   Find unique genotypes in sample based on total CN
            vec_unique_totCN_genotypes <- c(vec_unique_genotypes[1])
            if (length(vec_unique_genotypes) >= 2) {
                for (i in 1:length(vec_unique_genotypes)) {
                    genotype_next <- vec_unique_genotypes[i]
                    flag_new <- 1
                    for (j in 1:length(vec_unique_totCN_genotypes)) {
                        genotype_old <- vec_unique_totCN_genotypes[j]
                        flag_identical <- 1
                        if (identical(genotype_list_ploidy_chrom[[genotype_next]], genotype_list_ploidy_chrom[[genotype_old]]) == FALSE) {
                            flag_identical <- 0
                        } else {
                            for (chrom in 1:N_chromosomes) {
                                for (strand in 1:genotype_list_ploidy_chrom[[genotype_next]][chrom]) {
                                    if (identical(genotype_list_ploidy_block[[genotype_next]][[chrom]][[strand]], genotype_list_ploidy_block[[genotype_old]][[chrom]][[strand]]) == FALSE) {
                                        flag_identical <- 0
                                        break
                                    }
                                }
                                if (flag_identical == 0) {
                                    break
                                }
                            }
                        }
                        if (flag_identical == 1) {
                            flag_new <- 0
                            break
                        }
                    }
                    if (flag_new == 1) {
                        vec_unique_totCN_genotypes <- c(vec_unique_totCN_genotypes, genotype_next)
                    }
                }
            }
            stats_sample$total_cn_based_n_clones[sample] <- length(vec_unique_totCN_genotypes)
        }
        #---Statistics of CNA events found in sample(s)
        stats_sample$perc_initial <- 0
        stats_sample$perc_driver_mut <- 0
        stats_sample$perc_WGD <- 0
        stats_sample$perc_missegregation <- 0
        evolution_origin <- simulation$clonal_evolution$evolution_origin
        evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
        for (sample in 1:nrow(stats_sample)) {
            vec_loc <- which(simulation$sample$all_sample_ID == stats_sample$Sample_ID[sample])
            vec_unique_genotypes <- unique(simulation$sample$all_sample_genotype[vec_loc])
            vec_unique_genotypes_pop <- rep(0, length = length(vec_unique_genotypes))
            for (clone in 1:length(vec_unique_genotypes)) {
                vec_unique_genotypes_pop[clone] <- length(which(simulation$sample$all_sample_genotype[vec_loc] == vec_unique_genotypes[clone]))
            }
            #   Find list of events in each genotype
            vec_unique_genotypes_event_list <- vector("list", length = length(vec_unique_genotypes))
            for (clone in 1:length(vec_unique_genotypes)) {
                genotype <- vec_unique_genotypes[clone]
                genotype_event_list <- c()
                while (genotype > 0) {
                    event_list <- evolution_genotype_changes[[genotype]]
                    for (event in 1:length(event_list)) {
                        genotype_event_list <- c(genotype_event_list, event_list[[event]][1])
                    }
                    genotype <- evolution_origin[genotype]
                }
                if (length(genotype_event_list) > 0) {
                    vec_unique_genotypes_event_list[[clone]] <- genotype_event_list
                }
            }
            #---Find sample percentages for each event group: initial/WGD/missegregation/...
            N_initial <- 0
            N_driver_mut <- 0
            N_WGD <- 0
            N_missegregation <- 0
            for (clone in 1:length(vec_unique_genotypes)) {
                clone_event_list <- vec_unique_genotypes_event_list[[clone]]
                if (length(clone_event_list) == 0) {
                    N_initial <- N_initial + vec_unique_genotypes_pop[clone]
                }
                if ("new-driver" %in% clone_event_list) {
                    N_driver_mut <- N_driver_mut + vec_unique_genotypes_pop[clone]
                }
                if ("whole-genome-duplication" %in% clone_event_list) {
                    N_WGD <- N_WGD + vec_unique_genotypes_pop[clone]
                }
                if ("missegregation" %in% clone_event_list) {
                    N_missegregation <- N_missegregation + vec_unique_genotypes_pop[clone]
                }
            }




            #---Find mean count of missegregations before and after WGD (weighted by clone populations)
            vec_WGD_preceded_by_miss_count <- c()
            vec_WGD_followed_by_miss_count <- c()
            vec_WGD_pop <- c()
            #   Find counts of missegregations before and after WGD
            for (clone in 1:length(vec_unique_genotypes)) {
                clone_event_list <- vec_unique_genotypes_event_list[[clone]]
                #   Only perform on clones having WGD
                if ("whole-genome-duplication" %in% clone_event_list) {
                    loc <- which(clone_event_list == "whole-genome-duplication")[1]
                    #   Extract events before WGD
                    if (loc > 1) {
                        clone_event_list_before_WGD <- clone_event_list[1:(loc - 1)]
                    } else {
                        clone_event_list_before_WGD <- c()
                    }
                    #   Extract events after WGD
                    if (loc < length(clone_event_list)) {
                        clone_event_list_after_WGD <- clone_event_list[(loc + 1):length(clone_event_list)]
                    } else {
                        clone_event_list_after_WGD <- c()
                    }
                    #   Count missegregations before WGD
                    vec_WGD_preceded_by_miss_count <- c(vec_WGD_preceded_by_miss_count, length(which(clone_event_list_before_WGD == "missegregation")))
                    #   Count missegregations after WGD
                    vec_WGD_followed_by_miss_count <- c(vec_WGD_followed_by_miss_count, length(which(clone_event_list_after_WGD == "missegregation")))
                    #   Store population of this WGD clone
                    vec_WGD_pop <- c(vec_WGD_pop, vec_unique_genotypes_pop[clone])
                }
            }
            #   Find average counts of missegregations before and after WGD (weighted by clone populations)
            if (length(vec_WGD_pop) == 0) {
                mean_miss_before_WGD <- NA
                mean_miss_after_WGD <- NA
            } else {
                mean_miss_before_WGD <- sum(vec_WGD_preceded_by_miss_count * vec_WGD_pop) / sum(vec_WGD_pop)
                mean_miss_after_WGD <- sum(vec_WGD_followed_by_miss_count * vec_WGD_pop) / sum(vec_WGD_pop)
            }





            #   Store statistics of CNA events found in sample(s)
            stats_sample$perc_initial[sample] <- 100 * N_initial / sum(vec_unique_genotypes_pop)
            stats_sample$perc_driver_mut[sample] <- 100 * N_driver_mut / sum(vec_unique_genotypes_pop)
            stats_sample$perc_WGD[sample] <- 100 * N_WGD / sum(vec_unique_genotypes_pop)
            stats_sample$perc_missegregation[sample] <- 100 * N_missegregation / sum(vec_unique_genotypes_pop)



            stats_sample$mean_miss_before_WGD <- mean_miss_before_WGD
            stats_sample$mean_miss_after_WGD <- mean_miss_after_WGD
        }
        #---Save dataframe for statistics
        simulation$statistics$stats_sample <- stats_sample
    }
    return(simulation)
}
