simulator_full_program <- function(model = "",
                                   n_simulations = 0,
                                   stage_final = 0,
                                   n_clones_min = 0,
                                   n_clones_max = Inf,
                                   save_simulation = FALSE,
                                   internal_nodes_cn_info = FALSE,
                                   save_newick_tree = FALSE,
                                   save_cn_profile = FALSE,
                                   format_cn_profile = "both",
                                   model_readcount = FALSE,
                                   pseudo_corrected_readcount = FALSE,
                                   apply_HMM = FALSE,
                                   apply_UMAP_on_HMM = FALSE,
                                   report_progress = TRUE,
                                   compute_parallel = FALSE,
                                   seed = Inf) {
    # ==================================OVERRIDE PARAMETERS IF NECESSARY
    if (apply_HMM == TRUE) {
        save_simulation <- TRUE
    }
    # =================CREATE WORKSPACE DIRECTORY FOR CN INFERENCE WORKS
    if (model_readcount == TRUE) {
        dir.create(model)
    }
    # ======================================MAIN LOOP OF CANCERSIMULATOR
    if (seed == Inf) {
        set.seed(Sys.time())
    } else {
        set.seed(seed)
    }
    if (compute_parallel == FALSE) {
        #-------------------------Run CancerSimulator in sequential mode
        all_simulations <- vector("list", length = n_simulations)
        for (iteration in 1:n_simulations) {
            if (report_progress == TRUE) {
                cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
                cat(paste("BEGINNING SIMULATION-", iteration, "...\n", sep = ""))
            }
            simulation <- one_simulation(model, stage_final, save_cn_profile, internal_nodes_cn_info, format_cn_profile, model_readcount, report_progress)
            all_simulations[[iteration]] <- simulation
        }
    } else {
        #---------------------------Run CancerSimulator in parallel mode
        #   Start parallel cluster
        numCores <- detectCores()
        cl <- makePSOCKcluster(numCores - 1)
        setDefaultCluster(cl)
        #   Prepare input parameters for CancerSimulator
        model_name <<- model_name
        n_clones_min <<- n_clones_min
        n_clones_max <<- n_clones_max
        model <<- model
        stage_final <<- stage_final
        save_cn_profile <<- save_cn_profile
        internal_nodes_cn_info <<- internal_nodes_cn_info
        format_cn_profile <<- format_cn_profile
        model_readcount <<- model_readcount
        report_progress <<- report_progress
        clusterExport(NULL, varlist = c(
            "model", "stage_final", "save_cn_profile", "internal_nodes_cn_info", "format_cn_profile", "model_readcount", "report_progress", "model_name", "n_clones_min", "n_clones_max",
            "one_simulation",
            "SIMULATOR_VARIABLES_for_simulation",
            "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
            "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
            "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
            "SIMULATOR_FULL_PHASE_2_main", "SIMULATOR_FULL_PHASE_3_main",
            "get_cn_profile", "rbindlist", "p2_cn_profiles_long", "p2_reads_from_cn"
        ))
        #   Run CancerSimulator in parallel
        all_simulations <- parLapply(NULL, 1:n_simulations, function(i) {
            simulation <- one_simulation(model, stage_final, save_cn_profile, internal_nodes_cn_info, format_cn_profile, model_readcount, report_progress)
            return(simulation)
        })
        #   Stop parallel cluster
        stopCluster(cl)
    }
    # =====================================PROCESSING OF CANCERSIMULATOR
    for (iteration in 1:n_simulations) {
        simulation <- all_simulations[[iteration]]
        #------------------------------------Save the simulation package
        if (save_simulation == TRUE) {
            if (report_progress == TRUE) cat("\nSave simulation package...\n")
            save(simulation, file = paste(model, "_simulation_", iteration, ".rda", sep = ""))
        }
        #-----------------------------Save GC & mappability as WIG files
        if (model_readcount == TRUE) {
            if (report_progress == TRUE) cat("\nSave GC & mappability in WIG format...\n")
            p2_write_gc_map_as_wig(filename_gc = paste(model, "_gc.wig", sep = ""), filename_map = paste(model, "_map.wig", sep = ""))
        }
        #-----------------------------------Save the sampled CN profiles
        if (save_cn_profile == TRUE) {
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nSave true CN profiles in long format...\n")
                cn_profiles_long <- simulation$sample$cn_profiles_long
                filename <- paste(model, "_cn_profiles_long_", iteration, ".csv", sep = "")
                write.csv(cn_profiles_long, filename, row.names = FALSE)
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nSave true CN profiles in wide format...\n")
                cn_profiles_wide <- simulation$sample$cn_profiles_wide
                filename <- paste(model, "_cn_profiles_wide_", iteration, ".csv", sep = "")
                write.csv(cn_profiles_wide, filename, row.names = FALSE)
            }
        }
        #------------------Save the table of CN events in cell phylogeny
        if (internal_nodes_cn_info == TRUE) {
            if (report_progress == TRUE) cat("\nSave table of CN events in cell phylogeny...\n")
            hclust_CN_events <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_CN_events
            filename <- paste(model, "_cn_events_", iteration, ".csv", sep = "")
            write.csv(hclust_CN_events, filename, row.names = FALSE)
        }
        #-------------------------------------Save the noisy CN profiles
        if (model_readcount == TRUE) {
            if (report_progress == TRUE) cat("\nSave noisy CN profiles in long format...\n")
            noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
            filename <- paste(model, "_noisy_cn_profiles_long_", iteration, ".csv", sep = "")
            write.csv(noisy_cn_profiles_long, filename, row.names = FALSE)
            if (report_progress == TRUE) cat("\nSave noisy CN profiles in WIG format...\n")
            sample_cell_ID <- simulation$sample$sample_cell_ID
            noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
            for (cell in 1:length(sample_cell_ID)) {
                cell_ID <- sample_cell_ID[cell]
                filename <- paste(model, "/", model, "_noisy_cn_profiles_long_", iteration, "_", cell_ID, ".wig", sep = "")
                p2_write_cn_as_wig(filename, noisy_cn_profiles_long, cell_ID)
            }
        }
        #-------------------------Save the sampled cells' phylogeny tree
        if (save_newick_tree == TRUE) {
            if (report_progress == TRUE) cat("\nSave sampled cells' phylogeny in Newick format...\n")
            cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
            filename <- paste(model, "_cell_phylogeny_", iteration, ".newick", sep = "")
            write(hc2Newick_MODIFIED(cell_phylogeny_hclust), file = filename)
        }
        #----------------------------------Save statistics of simulation
        if (stage_final >= 1) {
            #   True clone counts in whole population over time
            if (iteration == 1) {
                true_n_clones_whole_pop_times <- simulation$statistics$true_n_clones_whole_pop_times
                true_n_clones_whole_pop_values <- matrix(simulation$statistics$true_n_clones_whole_pop_values, nrow = n_simulations, ncol = length(simulation$statistics$true_n_clones_whole_pop_values))
            } else {
                true_n_clones_whole_pop_values[iteration, ] <- simulation$statistics$true_n_clones_whole_pop_values
            }
        }
        if (stage_final >= 2) {
            #   True clone counts in sample
            if (iteration == 1) {
                true_n_clones_sample <- data.frame(matrix(simulation$statistics$stats_sample$true_n_clones, nrow = 1))
                colnames(true_n_clones_sample) <- simulation$statistics$stats_sample$Sample_ID
            } else {
                true_n_clones_sample[nrow(true_n_clones_sample) + 1, ] <- simulation$statistics$stats_sample$true_n_clones
            }
            #   Total-CN-based clone counts in sample
            if (iteration == 1) {
                total_cn_based_n_clones_sample <- data.frame(matrix(simulation$statistics$stats_sample$total_cn_based_n_clones, nrow = 1))
                colnames(total_cn_based_n_clones_sample) <- simulation$statistics$stats_sample$Sample_ID
            } else {
                total_cn_based_n_clones_sample[nrow(total_cn_based_n_clones_sample) + 1, ] <- simulation$statistics$stats_sample$total_cn_based_n_clones
            }
            #   Event class percentages in sample
            if (iteration == 1) {
                perc_initial_sample <- data.frame(matrix(simulation$statistics$stats_sample$perc_initial, nrow = 1))
                colnames(perc_initial_sample) <- simulation$statistics$stats_sample$Sample_ID
                perc_driver_mut_sample <- data.frame(matrix(simulation$statistics$stats_sample$perc_driver_mut, nrow = 1))
                colnames(perc_driver_mut_sample) <- simulation$statistics$stats_sample$Sample_ID
                perc_WGD_sample <- data.frame(matrix(simulation$statistics$stats_sample$perc_WGD, nrow = 1))
                colnames(perc_WGD_sample) <- simulation$statistics$stats_sample$Sample_ID
                perc_missegregation_sample <- data.frame(matrix(simulation$statistics$stats_sample$perc_missegregation, nrow = 1))
                colnames(perc_missegregation_sample) <- simulation$statistics$stats_sample$Sample_ID
            } else {
                perc_initial_sample[nrow(perc_initial_sample) + 1, ] <- simulation$statistics$stats_sample$perc_initial
                perc_driver_mut_sample[nrow(perc_driver_mut_sample) + 1, ] <- simulation$statistics$stats_sample$perc_driver_mut
                perc_WGD_sample[nrow(perc_WGD_sample) + 1, ] <- simulation$statistics$stats_sample$perc_WGD
                perc_missegregation_sample[nrow(perc_missegregation_sample) + 1, ] <- simulation$statistics$stats_sample$perc_missegregation
            }
        }
    }
    # ===============================================DOWNSTREAM ANALYSIS
    #----------------------------------Compute statistics of simulations
    simulation_statistics <- list()
    if (stage_final >= 1) {
        #   True clone counts in whole population
        if (standard_time_unit == "year") {
            true_n_clones_whole_pop_times <- true_n_clones_whole_pop_times / 365
        } else {
            if (standard_time_unit == "month") {
                true_n_clones_whole_pop_times <- true_n_clones_whole_pop_times / 30
            }
        }
        true_n_clones_whole_pop_mean <- colMeans(true_n_clones_whole_pop_values)
        true_n_clones_whole_pop_sd <- colSds(true_n_clones_whole_pop_values)
        simulation_statistics$true_n_clones_whole_pop_times <- true_n_clones_whole_pop_times
        simulation_statistics$true_n_clones_whole_pop_values <- true_n_clones_whole_pop_values
        simulation_statistics$true_n_clones_whole_pop_mean <- true_n_clones_whole_pop_mean
        simulation_statistics$true_n_clones_whole_pop_sd <- true_n_clones_whole_pop_sd
    }
    if (stage_final >= 2) {
        statistics_sample <- data.frame(matrix(0, nrow = 12, ncol = ncol(true_n_clones_sample)))
        colnames(statistics_sample) <- colnames(true_n_clones_sample)
        rownames(statistics_sample) <- c(
            "true_n_clones_mean", "true_n_clones_sd",
            "total_cn_based_n_clones_mean", "total_cn_based_n_clones_sd",
            "perc_initial_mean", "perc_initial_sd",
            "perc_driver_mut_mean", "perc_driver_mut_sd",
            "perc_WGD_mean", "perc_WGD_sd",
            "perc_missegregation_mean", "perc_missegregation_sd"
        )
        for (sample in 1:ncol(true_n_clones_sample)) {
            #   True clone counts in sample
            statistics_sample[1, sample] <- mean(true_n_clones_sample[, sample])
            statistics_sample[2, sample] <- sd(true_n_clones_sample[, sample])
            #   Total-CN-based clone counts in sample
            statistics_sample[3, sample] <- mean(total_cn_based_n_clones_sample[, sample])
            statistics_sample[4, sample] <- sd(total_cn_based_n_clones_sample[, sample])
            #   Percentages of sample in each event class
            statistics_sample[5, sample] <- mean(perc_initial_sample[, sample])
            statistics_sample[6, sample] <- sd(perc_initial_sample[, sample])
            statistics_sample[7, sample] <- mean(perc_driver_mut_sample[, sample])
            statistics_sample[8, sample] <- sd(perc_driver_mut_sample[, sample])
            statistics_sample[9, sample] <- mean(perc_WGD_sample[, sample])
            statistics_sample[10, sample] <- sd(perc_WGD_sample[, sample])
            statistics_sample[11, sample] <- mean(perc_missegregation_sample[, sample])
            statistics_sample[12, sample] <- sd(perc_missegregation_sample[, sample])
        }
        simulation_statistics$true_n_clones_sample <- true_n_clones_sample
        simulation_statistics$total_cn_based_n_clones_sample <- total_cn_based_n_clones_sample
        simulation_statistics$perc_initial_sample <- perc_initial_sample
        simulation_statistics$perc_driver_mut_sample <- perc_driver_mut_sample
        simulation_statistics$perc_WGD_sample <- perc_WGD_sample
        simulation_statistics$perc_missegregation_sample <- perc_missegregation_sample
        simulation_statistics$statistics_sample <- statistics_sample
    }
    #-----------------------Infer CN from noisy readcounts using HMMcopy
    if (apply_HMM == TRUE) {
        if (report_progress == TRUE) {
            cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
            cat("RUNNING HMMCOPY FOR ALL SIMULATIONS\n")
        }
        #   Run HMMcopy for each individual cell
        if (compute_parallel == TRUE) {
            flag_parallel <- 1
        } else {
            flag_parallel <- 0
        }
        system("chmod -x hmmcopy.bash")
        system(paste("bash hmmcopy.bash ", model, " ", n_simulations, " ", flag_parallel, sep = ""))
        #   Append the simulation package RDA with HMMcopy inference
        append_with_hmm(model = model, n_simulations = n_simulations, UMAP = apply_UMAP_on_HMM, pseudo_corrected_readcount = pseudo_corrected_readcount)
    }
    # =================================================OUTPUT SIMULATION
    return(simulation_statistics)
}

one_simulation <- function(model, stage_final, save_cn_profile, internal_nodes_cn_info, format_cn_profile, model_readcount, report_progress) {
    # =============================================LOAD MODEL PARAMETERS
    SIMULATOR_VARIABLES_for_simulation(model)
    # ============================================PRODUCE ONE SIMULATION
    flag_success <- 0
    while (flag_success == 0) {
        #----------------------------------Simulate the clonal evolution
        if (report_progress == TRUE) cat("\nStage 1: clonal evolution...\n")
        output <- SIMULATOR_FULL_PHASE_1_main(report_progress)
        flag_success <- output$flag_success
        package_clonal_evolution <- output$package_clonal_evolution
        if (flag_success == 0) {
            cat("SIMULATION CONDITIONS NOT SATISFIED; REDOING...\n")
            next
        }
        #--------------------------------------------Simulate the sample
        if (stage_final >= 2) {
            if (report_progress == TRUE) cat("\nStage 2: sampling...\n")
            package_sample <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution, report_progress)
            N_clones <- nrow(package_sample$table_clone_ID_vs_letters)
            if ((N_clones < n_clones_min) | (N_clones > n_clones_max)) {
                flag_success <- 0
                cat("SIMULATION CONDITIONS NOT SATISFIED; REDOING...\n")
                next
            }
        }
        #---------------------------Simulate the phylogeny of the sample
        if (stage_final >= 3) {
            if (report_progress == TRUE) cat("\nStage 3: sample phylogeny...\n")
            package_sample_phylogeny <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution, package_sample)
        }
    }
    # =============================COMPUTE STATISTICS FOR THE SIMULATION



    # ======================PREPARE DATA FROM PHASE 1 (CLONAL EVOLUTION)
    if (stage_final >= 1) {
        simulation <- list()
        simulation$clonal_evolution <- package_clonal_evolution
        simulation$statistics <- list()
        #---Compute statistics - true clone counts in whole population
        true_n_clones_whole_pop_times <- simulation$clonal_evolution$evolution_traj_time
        true_n_clones_whole_pop_values <- rep(0, length = length(true_n_clones_whole_pop_times))
        for (step in 1:length(true_n_clones_whole_pop_values)) {
            true_n_clones_whole_pop_values[step] <- length(simulation$clonal_evolution$evolution_traj_clonal_ID[[step]])
        }
        simulation$statistics$true_n_clones_whole_pop_times <- true_n_clones_whole_pop_times
        simulation$statistics$true_n_clones_whole_pop_values <- true_n_clones_whole_pop_values
    }
    # ==============================PREPARE DATA FROM PHASE 2 (SAMPLING)
    if (stage_final >= 2) {
        simulation$sample <- package_sample
        #---Build CN profile tables in long or wide format or both
        if (save_cn_profile == TRUE) {
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nExtra: build CN profile table in long format...\n")
                simulation <- p2_cn_profiles_long(simulation)
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nExtra: build CN profile table in wide format...\n")
                simulation <- p2_cn_profiles_wide(simulation)
            }
        } else {
            if (model_readcount == TRUE) {
                if (report_progress == TRUE) cat("\nExtra: build CN profile table in long format (for modeling CN with noise)...\n")
                simulation <- p2_cn_profiles_long(simulation)
            }
        }
        #---Simulate noisy readcounts
        if (model_readcount == TRUE) {
            if (report_progress == TRUE) cat("\nExtra: simulate CN profiles with noise...\n")
            simulation <- p2_reads_from_cn(simulation, report_progress)
        }
        #---Compute statistics - prepare dataframe for statistics
        stats_sample <- Table_sampling
        #---Compute statistics - true clone counts in sample
        stats_sample$true_n_clones <- 0
        for (sample in 1:nrow(stats_sample)) {
            vec_loc <- which(simulation$sample$all_sample_ID == stats_sample$Sample_ID[sample])
            stats_sample$true_n_clones[sample] <- length(unique(simulation$sample$all_sample_genotype[vec_loc]))
        }
        #---Compute statistics - clone counts in sample based on total CN
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
        #---Compute statistics - sample percentages that are initial/WGD/missegregation/...
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
                        genotype_event_list <- unique(c(genotype_event_list, event_list[[event]][1]))
                    }
                    genotype <- evolution_origin[genotype]
                }
                if (length(genotype_event_list) > 0) {
                    vec_unique_genotypes_event_list[[clone]] <- genotype_event_list
                }
            }
            #   Find sample percentages for each event group
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
            stats_sample$perc_initial[sample] <- 100 * N_initial / sum(vec_unique_genotypes_pop)
            stats_sample$perc_driver_mut[sample] <- 100 * N_driver_mut / sum(vec_unique_genotypes_pop)
            stats_sample$perc_WGD[sample] <- 100 * N_WGD / sum(vec_unique_genotypes_pop)
            stats_sample$perc_missegregation[sample] <- 100 * N_missegregation / sum(vec_unique_genotypes_pop)
        }
        #---Compute statistics - save dataframe for statistics
        simulation$statistics$stats_sample <- stats_sample
    }
    # ======================PREPARE DATA FROM PHASE 3 (SAMPLE PHYLOGENY)
    if (stage_final >= 3) {
        simulation$sample_phylogeny <- package_sample_phylogeny
        #---Supplement sample phylogeny data with internal nodes
        if (save_cn_profile == TRUE) {
            if (internal_nodes_cn_info == TRUE) {
                if (report_progress == TRUE) cat("\nExtra: find CN information for internal nodes...\n")
                simulation <- p3_cn_profiles_internal(simulation)
                if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                    if (report_progress == TRUE) cat("\nExtra: add CN profiles for internal nodes in long format...\n")
                    simulation <- p3_internal_node_cn_profiles_long(simulation)
                }
                if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                    if (report_progress == TRUE) cat("\nExtra: build CN profile table in wide format...\n")
                    simulation <- p3_internal_node_cn_profiles_wide(simulation)
                }
                if (report_progress == TRUE) cat("\nExtra: build table of CN events...\n")
                simulation <- p3_cn_events_table(simulation)
            }
        }
    }
    # =====================================OUTPUT THE SIMULATION PACKAGE
    return(simulation)
}
