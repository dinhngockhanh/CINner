#' Produce simulations and outputs files as requested by user
#'
#' @param model ...
#' @param n_simulations ...
#' @param stage_final ...
#' @param n_clones_min ...
#' @param n_clones_max ...
#' @param save_simulation ...
#' @param internal_nodes_cn_info ...
#' @param save_newick_tree ...
#' @param save_cn_profile ...
#' @param save_cn_clones ...
#' @param format_cn_profile ...
#' @param model_readcount ...
#' @param pseudo_corrected_readcount ...
#' @param apply_HMM ...
#' @param apply_UMAP_on_HMM ...
#' @param report_progress ...
#' @param compute_parallel ...
#' @param seed ...
#' @export
simulator_full_program <- function(model = "",
                                   model_prefix = "",
                                   n_simulations = 0,
                                   stage_final = 0,
                                   n_clones_min = 0,
                                   n_clones_max = Inf,
                                   save_simulation = TRUE,
                                   internal_nodes_cn_info = FALSE,
                                   save_newick_tree = FALSE,
                                   save_cn_profile = FALSE,
                                   save_cn_clones = FALSE,
                                   format_cn_profile = "both",
                                   model_readcount = FALSE,
                                   pseudo_corrected_readcount = FALSE,
                                   apply_HMM = FALSE,
                                   apply_UMAP_on_HMM = FALSE,
                                   report_progress = TRUE,
                                   compute_parallel = FALSE,
                                   seed = Inf,
                                   output_variables = c(),
                                   n_cores = NULL) {
    if (class(model) == "character") {
        model_parameters <- model
        model_prefix <- model
    } else if (class(model) == "list") {
        model_parameters <- model
    }
    # ==================================OVERRIDE PARAMETERS IF NECESSARY
    if (apply_HMM == TRUE) {
        save_simulation <- TRUE
    }
    # =================CREATE WORKSPACE DIRECTORY FOR CN INFERENCE WORKS
    if (model_readcount == TRUE) {
        dir.create(model_prefix)
    }
    # ======================================MAIN LOOP OF CANCERSIMULATOR
    if (seed == Inf) {
        set.seed(Sys.time())
    } else {
        set.seed(seed)
    }
    if (compute_parallel == FALSE) {
        #-------------------------Run CancerSimulator in sequential mode
        many_sims <- list()
        for (iteration in 1:n_simulations) {
            if (report_progress == TRUE) {
                cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
                cat(paste("BEGINNING SIMULATION-", iteration, "...\n", sep = ""))
            }
            many_sims[[iteration]] <- one_simulation(
                iteration,
                model_parameters,
                model_prefix,
                stage_final,
                n_clones_min,
                n_clones_max,
                save_simulation,
                internal_nodes_cn_info,
                save_newick_tree,
                save_cn_profile,
                save_cn_clones,
                format_cn_profile,
                model_readcount,
                report_progress,
                output_variables
            )
        }
    } else {
        #---------------------------Run CancerSimulator in parallel mode
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        # setDefaultCluster(cl)
        if (report_progress == TRUE) {
            cat(paste("\nSTARTED PARALLEL CLUSTER WITH ", numCores - 1, " CORES...\n", sep = ""))
        }
        #   Prepare input parameters for CancerSimulator
        model_parameters <<- model_parameters
        model_prefix <<- model_prefix
        stage_final <<- stage_final
        n_clones_min <<- n_clones_min
        n_clones_max <<- n_clones_max
        save_simulation <<- save_simulation
        internal_nodes_cn_info <<- internal_nodes_cn_info
        save_newick_tree <<- save_newick_tree
        save_cn_profile <<- save_cn_profile
        save_cn_clones <<- save_cn_clones
        format_cn_profile <<- format_cn_profile
        model_readcount <<- model_readcount
        report_progress <<- report_progress
        output_variables <<- output_variables
        clusterExport(cl, varlist = c(
            "model_parameters",
            "model_prefix",
            "stage_final",
            "n_clones_min",
            "n_clones_max",
            "save_simulation",
            "internal_nodes_cn_info",
            "save_newick_tree",
            "save_cn_profile",
            "save_cn_clones",
            "format_cn_profile",
            "model_readcount",
            "report_progress",
            "one_simulation",
            "output_variables",
            "hc2Newick_MODIFIED", "hc2Newick",
            "SIMULATOR_VARIABLES_for_simulation",
            "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
            "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
            "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
            "SIMULATOR_FULL_PHASE_2_main", "SIMULATOR_FULL_PHASE_3_main",
            "get_cn_profile", "rbindlist", "p2_cn_profiles_long", "p2_readcount_model"
        ))
        library(ape)
        clusterEvalQ(cl = cl, require(ape))
        #   Run CancerSimulator in parallel
        if (report_progress == TRUE) {
            many_sims <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
                one_sim <- one_simulation(
                    iteration,
                    model_parameters,
                    model_prefix,
                    stage_final,
                    n_clones_min,
                    n_clones_max,
                    save_simulation,
                    internal_nodes_cn_info,
                    save_newick_tree,
                    save_cn_profile,
                    save_cn_clones,
                    format_cn_profile,
                    model_readcount,
                    report_progress,
                    output_variables
                )
                return(one_sim)
            })
        } else {
            many_sims <- parLapply(cl, 1:n_simulations, function(iteration) {
                one_sim <- one_simulation(
                    iteration,
                    model_parameters,
                    model_prefix,
                    stage_final,
                    n_clones_min,
                    n_clones_max,
                    save_simulation,
                    internal_nodes_cn_info,
                    save_newick_tree,
                    save_cn_profile,
                    save_cn_clones,
                    format_cn_profile,
                    model_readcount,
                    report_progress,
                    output_variables
                )
                return(one_sim)
            })
        }
        #   Stop parallel cluster
        stopCluster(cl)
    }
    # ================================SAVE GC & MAPPABILITY AS WIG FILES
    if (model_readcount == TRUE) {
        if (report_progress == TRUE) cat("\nSave GC & mappability in WIG format...\n")
        p2_write_gc_map_as_wig(filename_gc = paste(model_prefix, "_gc.wig", sep = ""), filename_map = paste(model_prefix, "_map.wig", sep = ""))
    }
    # ======================INFER CN FROM NOISY READCOUNTS USING HMMCOPY
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
        system(paste("bash hmmcopy.bash ", model_prefix, " ", n_simulations, " ", flag_parallel, sep = ""))
        #   Append the simulation package RDA with HMMcopy inference
        append_with_hmm(model = model, n_simulations = n_simulations, UMAP = apply_UMAP_on_HMM, pseudo_corrected_readcount = pseudo_corrected_readcount)
    }
    # ==================================OUTPUT COLLECTION OF SIMULATIONS
    return(many_sims)
}

#' @export
one_simulation <- function(iteration,
                           model_parameters,
                           model_prefix,
                           stage_final,
                           n_clones_min,
                           n_clones_max,
                           save_simulation,
                           internal_nodes_cn_info,
                           save_newick_tree,
                           save_cn_profile,
                           save_cn_clones,
                           format_cn_profile,
                           model_readcount,
                           report_progress,
                           output_variables) {
    # =============================================LOAD MODEL PARAMETERS
    SIMULATOR_VARIABLES_for_simulation(model_parameters)
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
            package_sample_phylogeny <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution, package_sample, report_progress)
        }
    }
    # ======================PREPARE DATA FROM PHASE 1 (CLONAL EVOLUTION)
    simulation <- list()
    if (stage_final >= 1) {
        simulation$clonal_evolution <- package_clonal_evolution
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
                if (report_progress == TRUE) cat("\nExtra: build CN profile table in long format (for modeling readcounts with noise & bias)...\n")
                simulation <- p2_cn_profiles_long(simulation)
            }
        }
        #---Simulate noisy readcounts
        if (model_readcount == TRUE) {
            if (report_progress == TRUE) cat("\nExtra: simulate readcount profiles with noise & bias...\n")
            simulation <- p2_readcount_model(simulation, report_progress)
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nExtra: build noisy & biased readcount profile table in wide format...\n")
                simulation <- p2_readcount_model_wide(simulation, report_progress)
            }
        }
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
    # ======================================OUTPUT FILES FROM SIMULATION
    #----------------------------------------Save the simulation package
    if (save_simulation == TRUE) {
        if (report_progress == TRUE) cat("\nSave simulation package...\n")
        save(simulation, file = paste(model_prefix, "_simulation_", iteration, ".rda", sep = ""))
    }
    #---------------------------------------Save the sampled CN profiles
    if (save_cn_profile == TRUE) {
        if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("\nSave true CN profiles in long format...\n")
            cn_profiles_long <- simulation$sample$cn_profiles_long
            filename <- paste(model_prefix, "_cn_profiles_long_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_long, filename, row.names = FALSE)
        }
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("\nSave true CN profiles in wide format...\n")
            cn_profiles_wide <- simulation$sample$cn_profiles_wide
            filename <- paste(model_prefix, "_cn_profiles_wide_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_wide, filename, row.names = FALSE)
        }
    }
    #---------------------------Save the clonal identities of every cell
    if (save_cn_clones == TRUE) {
        if (report_progress == TRUE) cat("\nSave table of cell-clone mapping...\n")
        table_cell_clone <- simulation$sample$table_cell_clone
        filename <- paste(model_prefix, "_cn_profiles_clonal_mapping_", iteration, ".csv", sep = "")
        write.csv(table_cell_clone, filename, row.names = FALSE)
    }
    #----------------------Save the table of CN events in cell phylogeny
    if (internal_nodes_cn_info == TRUE) {
        if (report_progress == TRUE) cat("\nSave table of CN events in cell phylogeny...\n")
        hclust_CN_events <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_CN_events
        filename <- paste(model_prefix, "_cn_events_", iteration, ".csv", sep = "")
        write.csv(hclust_CN_events, filename, row.names = FALSE)
    }
    #-----------------------------------------Save the noisy CN profiles
    if (model_readcount == TRUE) {
        if (save_cn_profile == TRUE) {
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nSave noisy CN profiles in long format...\n")
                noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
                filename <- paste(model_prefix, "_noisy_cn_profiles_long_", iteration, ".csv", sep = "")
                write.csv(noisy_cn_profiles_long, filename, row.names = FALSE)
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("\nSave noisy CN profiles in wide format...\n")
                noisy_cn_profiles_wide <- simulation$sample$noisy_cn_profiles_wide
                filename <- paste(model_prefix, "_noisy_cn_profiles_wide_", iteration, ".csv", sep = "")
                write.csv(noisy_cn_profiles_wide, filename, row.names = FALSE)
            }
        }
        if (report_progress == TRUE) cat("\nSave noisy CN profiles in WIG format...\n")
        sample_cell_ID <- simulation$sample$sample_cell_ID
        noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
        for (cell in 1:length(sample_cell_ID)) {
            cell_ID <- sample_cell_ID[cell]
            filename <- paste(model_prefix, "/", model_prefix, "_noisy_cn_profiles_long_", iteration, "_", cell_ID, ".wig", sep = "")
            p2_write_cn_as_wig(filename, noisy_cn_profiles_long, cell_ID)
        }
    }
    #-----------------------------Save the sampled cells' phylogeny tree
    if (save_newick_tree == TRUE) {
        if (report_progress == TRUE) cat("\nSave sampled cells' phylogeny in Newick format...\n")
        cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
        filename <- paste(model_prefix, "_cell_phylogeny_", iteration, ".newick", sep = "")
        write(hc2Newick_MODIFIED(cell_phylogeny_hclust), file = filename)

        clone_phylogeny_hclust <- simulation$sample_phylogeny$clone_phylogeny_hclust
        filename <- paste(model_prefix, "_clone_phylogeny_", iteration, ".newick", sep = "")
        write(hc2Newick(clone_phylogeny_hclust), file = filename)
        # write(hc2Newick_MODIFIED(clone_phylogeny_hclust), file = filename)
    }
    #---------------------------------------Return the simulation result
    if (length(output_variables) == 0) {
        simulation_output <- simulation
    } else {
        simulation_output <- list()
        if ("all_sample_genotype" %in% output_variables) {
            simulation_output$sample$all_sample_genotype <- simulation$sample$all_sample_genotype
        }
        if ("sample_cell_ID" %in% output_variables) {
            simulation_output$sample$sample_cell_ID <- simulation$sample$sample_cell_ID
        }
        if ("sample_genotype_unique" %in% output_variables) {
            simulation_output$sample$sample_genotype_unique <- simulation$sample$sample_genotype_unique
        }
        if ("sample_genotype_unique_profile" %in% output_variables) {
            simulation_output$sample$sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
        }
        if ("cell_phylogeny_hclust" %in% output_variables) {
            simulation_output$sample_phylogeny$cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
        }
    }
    return(simulation_output)
}
