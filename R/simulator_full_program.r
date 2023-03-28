#' Produce simulations and outputs files as requested by user
#'
#' @param model_tmp ...
#' @param n_simulations ...
#' @param stage_final ...
#' @param n_clones_min ...
#' @param n_clones_max ...
#' @param save_simulation ...
#' @param neutral_variations ...
#' @param internal_nodes_cn_info ...
#' @param save_newick_tree ...
#' @param save_cn_profile ...
#' @param save_cn_clones ...
#' @param format_cn_profile ...
#' @param model_readcount ...
#' @param pseudo_corrected_readcount ...
#' @param HMM ...
#' @param report_progress ...
#' @param compute_parallel ...
#' @param seed ...
#' @export
simulator_full_program <- function(model = "",
                                   model_prefix = "",
                                   n_simulations = 0,
                                   stage_final = 1,
                                   n_clones_min = 1,
                                   n_clones_max = Inf,
                                   save_simulation = TRUE,
                                   neutral_variations = FALSE,
                                   internal_nodes_cn_info = FALSE,
                                   save_newick_tree = FALSE,
                                   save_cn_profile = FALSE,
                                   save_cn_clones = FALSE,
                                   build_cn = TRUE,
                                   format_cn_profile = "long",
                                   model_readcount = FALSE,
                                   model_readcount_base = "all",
                                   pseudo_corrected_readcount = FALSE,
                                   HMM = FALSE,
                                   HMM_containner = "docker",
                                   folder_workplace = NULL,
                                   report_progress = TRUE,
                                   compute_parallel = FALSE,
                                   seed = Inf,
                                   output_variables = c(),
                                   n_cores = NULL,
                                   R_libPaths = NULL) {
    library(data.table)
    library(ape)
    library(ctc)

    model_tmp <- model
    model_prefix_tmp <- model_prefix
    folder_workplace_tmp <- folder_workplace

    if (class(model_tmp) == "character") {
        model_parameters <- model_tmp
        model_prefix_tmp <- model_tmp
    } else if (class(model_tmp) == "list") {
        model_parameters <- model_tmp
    }
    if (is.null(folder_workplace_tmp)) {
        folder_workplace_tmp <- ""
    } else {
        dir.create(folder_workplace_tmp)
        folder_workplace_tmp <- paste(folder_workplace_tmp, "/", sep = "")
    }
    # ==================================OVERRIDE PARAMETERS IF NECESSARY
    if (pseudo_corrected_readcount == TRUE) {
        stage_final <- max(stage_final, 2)
        HMM <- TRUE
    }
    if (HMM == TRUE) {
        stage_final <- max(stage_final, 2)
        model_readcount <- TRUE
        save_simulation <- TRUE
    }
    if (save_cn_profile == TRUE) stage_final <- max(stage_final, 2)
    if (save_cn_clones == TRUE) stage_final <- max(stage_final, 2)
    if (model_readcount == TRUE) stage_final <- max(stage_final, 2)

    if (internal_nodes_cn_info == TRUE) stage_final <- max(stage_final, 3)
    if (save_newick_tree == TRUE) stage_final <- max(stage_final, 3)

    if (neutral_variations == TRUE) stage_final <- max(stage_final, 4)
    # ============PREPARE WORKSPACE DIRECTORY AND IMPORT NECESSARY FILES
    if (HMM == TRUE) {
        dir.create(model_prefix_tmp)
        if (HMM_containner == "singularity") {
            if (file.exists("hmmcopy_v0.0.45.sif") == FALSE) {
                cat("INSTALLING SINGULARITY FOR HMMCOPY-V0.0.45...\n")
                system("singularity pull docker://quay.io/mondrianscwgs/hmmcopy:v0.0.45")
            }
            file.copy("hmmcopy_v0.0.45.sif", paste(model_prefix_tmp, "/hmmcopy_v0.0.45.sif", sep = ""))
        } else if (HMM_containner == "docker") {
            system("docker pull quay.io/mondrianscwgs/hmmcopy:v0.0.45")
        }
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
                model_prefix_tmp,
                stage_final,
                n_clones_min,
                n_clones_max,
                save_simulation,
                build_cn,
                neutral_variations,
                internal_nodes_cn_info,
                save_newick_tree,
                save_cn_profile,
                save_cn_clones,
                format_cn_profile,
                model_readcount,
                model_readcount_base,
                pseudo_corrected_readcount,
                HMM,
                HMM_containner,
                folder_workplace_tmp,
                report_progress,
                output_variables
            )
        }
    } else {
        library(parallel)
        library(ctc)
        library(ape)
        #---------------------------Run CancerSimulator in parallel mode
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        #   Configure parallel pool
        cl <- makePSOCKcluster(numCores - 1)
        if (report_progress == TRUE) {
            cat(paste("\nParallel cluster with ", numCores - 1, " cores...\n", sep = ""))
        }
        if (is.null(R_libPaths) == FALSE) {
            R_libPaths <<- R_libPaths
            clusterExport(cl, varlist = c("R_libPaths"))
            clusterEvalQ(cl = cl, .libPaths(R_libPaths))
        }
        clusterEvalQ(cl = cl, library(ape))
        clusterEvalQ(cl = cl, library(ctc))
        #   Prepare input parameters for CancerSimulator
        model_parameters <<- model_parameters
        model_prefix_tmp <<- model_prefix_tmp
        stage_final <<- stage_final
        n_clones_min <<- n_clones_min
        n_clones_max <<- n_clones_max
        save_simulation <<- save_simulation
        build_cn <<- build_cn
        neutral_variations <<- neutral_variations
        internal_nodes_cn_info <<- internal_nodes_cn_info
        save_newick_tree <<- save_newick_tree
        save_cn_profile <<- save_cn_profile
        save_cn_clones <<- save_cn_clones
        format_cn_profile <<- format_cn_profile
        model_readcount <<- model_readcount
        model_readcount_base <<- model_readcount_base
        pseudo_corrected_readcount <<- pseudo_corrected_readcount
        HMM <<- HMM
        HMM_containner <<- HMM_containner
        folder_workplace_tmp <<- folder_workplace_tmp
        report_progress <<- report_progress
        output_variables <<- output_variables
        clusterExport(cl, varlist = c(
            "model_parameters",
            "model_prefix_tmp",
            "stage_final",
            "n_clones_min",
            "n_clones_max",
            "save_simulation",
            "build_cn",
            "neutral_variations",
            "internal_nodes_cn_info",
            "save_newick_tree",
            "save_cn_profile",
            "save_cn_clones",
            "format_cn_profile",
            "model_readcount",
            "model_readcount_base",
            "pseudo_corrected_readcount",
            "HMM",
            "HMM_containner",
            "folder_workplace_tmp",
            "report_progress",
            "one_simulation",
            "output_variables",
            "hc2Newick_MODIFIED", "hc2Newick",
            "SIMULATOR_VARIABLES_for_simulation",
            "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
            "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
            "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
            "SIMULATOR_FULL_PHASE_2_main", "SIMULATOR_FULL_PHASE_3_main", "SIMULATOR_FULL_PHASE_4_main",
            "get_cn_profile", "rbindlist",
            "p0_write_cn_as_wig", "p0_write_gc_map_as_wig", "p0_append_with_hmm",
            "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model", "p2_readcount_model_wide",
            "p3_cn_events_table", "p3_cn_profiles_internal", "p3_internal_node_cn_profiles_long", "p3_internal_node_cn_profiles_wide",
            "p4_cn_profiles_long", "p4_cn_profiles_wide", "p4_readcount_model", "p4_readcount_model_wide"
        ))
        #   Run CancerSimulator in parallel
        if (report_progress == TRUE) {
            library(pbapply)
            pbo <- pboptions(type = "txt")
            many_sims <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
                one_sim <- one_simulation(
                    iteration,
                    model_parameters,
                    model_prefix_tmp,
                    stage_final,
                    n_clones_min,
                    n_clones_max,
                    save_simulation,
                    build_cn,
                    neutral_variations,
                    internal_nodes_cn_info,
                    save_newick_tree,
                    save_cn_profile,
                    save_cn_clones,
                    format_cn_profile,
                    model_readcount,
                    model_readcount_base,
                    pseudo_corrected_readcount,
                    HMM,
                    HMM_containner,
                    folder_workplace_tmp,
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
                    model_prefix_tmp,
                    stage_final,
                    n_clones_min,
                    n_clones_max,
                    save_simulation,
                    build_cn,
                    neutral_variations,
                    internal_nodes_cn_info,
                    save_newick_tree,
                    save_cn_profile,
                    save_cn_clones,
                    format_cn_profile,
                    model_readcount,
                    model_readcount_base,
                    pseudo_corrected_readcount,
                    HMM,
                    HMM_containner,
                    folder_workplace_tmp,
                    report_progress,
                    output_variables
                )
                return(one_sim)
            })
        }
        #   Stop parallel cluster
        stopCluster(cl)
    }
    # ==================================OUTPUT COLLECTION OF SIMULATIONS
    return(many_sims)
}

#' @export
one_simulation <- function(iteration,
                           model_parameters,
                           model_prefix_tmp,
                           stage_final,
                           n_clones_min,
                           n_clones_max,
                           save_simulation,
                           build_cn,
                           neutral_variations,
                           internal_nodes_cn_info,
                           save_newick_tree,
                           save_cn_profile,
                           save_cn_clones,
                           format_cn_profile,
                           model_readcount,
                           model_readcount_base,
                           pseudo_corrected_readcount,
                           HMM,
                           HMM_containner,
                           folder_workplace_tmp,
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
        #------------------Simulate the neutral variations within clones
        if (stage_final >= 4) {
            if (report_progress == TRUE) cat("\nStage 4: subclonal neutral variations...\n")
            package_sample_with_neutral_variations <- SIMULATOR_FULL_PHASE_4_main(package_clonal_evolution, package_sample, package_sample_phylogeny, report_progress)
        }
    }
    if (report_progress == TRUE) cat("\n")
    # ======================PREPARE DATA FROM PHASE 1 (CLONAL EVOLUTION)
    simulation <- list()
    if (stage_final >= 1) {
        simulation$clonal_evolution <- package_clonal_evolution
    }
    # ==============================PREPARE DATA FROM PHASE 2 (SAMPLING)
    if (stage_final >= 2) {
        simulation$sample <- package_sample
        #---Build CN profile tables in long format
        if (build_cn == TRUE) {
            if (report_progress == TRUE) cat("Extra: build CN profile table in long format...\n")
            simulation <- p2_cn_profiles_long(simulation)
        }
        #---Build CN profile tables in wide format
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Extra: build CN profile table in wide format...\n")
            simulation <- p2_cn_profiles_wide(simulation)
        }
        #---Simulate noisy readcounts
        if ((model_readcount == TRUE) & ((model_readcount_base == "truth") | (model_readcount_base == "all"))) {
            if (report_progress == TRUE) cat("Extra: simulate readcount profiles with noise & bias in long format...\n")
            simulation <- p2_readcount_model(simulation, report_progress)
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("Extra: simulate readcount profiles with noise & bias in wide format...\n")
                simulation <- p2_readcount_model_wide(simulation, report_progress)
            }
        }
    }
    # ======================PREPARE DATA FROM PHASE 3 (SAMPLE PHYLOGENY)
    if (stage_final >= 3) {
        simulation$sample_phylogeny <- package_sample_phylogeny
        #---Supplement sample phylogeny data with internal nodes
        if (save_cn_profile == TRUE & internal_nodes_cn_info == TRUE) {
            if (report_progress == TRUE) cat("Extra: find CN information for internal nodes...\n")
            simulation <- p3_cn_profiles_internal(simulation)
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("Extra: add CN profiles for internal nodes in long format...\n")
                simulation <- p3_internal_node_cn_profiles_long(simulation)
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("Extra: build CN profile table in wide format...\n")
                simulation <- p3_internal_node_cn_profiles_wide(simulation)
            }
            if (report_progress == TRUE) cat("Extra: build table of CN events...\n")
            simulation <- p3_cn_events_table(simulation)
        }
    }
    # ===========================PREPARE DATA FROM PHASE 4 (NEUTRAL CNA)
    if (stage_final >= 4) {
        simulation$neutral_variations <- package_sample_with_neutral_variations
        #---Build CN profile tables with neutral variations in long format
        if (report_progress == TRUE) cat("Extra: build CN profile table with neutral variations in long format...\n")
        simulation <- p4_cn_profiles_long(simulation)
        #---Build CN profile tables with neutral variations in wide format
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Extra: build CN profile table with neutral variations in wide format...\n")
            simulation <- p4_cn_profiles_wide(simulation)
        }
        #---Simulate noisy readcounts
        if ((model_readcount == TRUE) & ((model_readcount_base == "neuvar") | (model_readcount_base == "all"))) {
            if (report_progress == TRUE) cat("Extra: simulate readcount profiles with neutral variations with noise & bias in long format...\n")
            simulation <- p4_readcount_model(simulation, report_progress)
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                if (report_progress == TRUE) cat("Extra: simulate readcount profiles with neutral variations with noise & bias in wide format...\n")
                simulation <- p4_readcount_model_wide(simulation, report_progress)
            }
        }
    }
    # ==================================PREPARE DATA FROM VARIOUS PHASES
    #-------------Save noisy CN profiles on base of true profiles in WIG
    if ((HMM == TRUE) & ((model_readcount_base == "truth") | (model_readcount_base == "all"))) {
        if (report_progress == TRUE) cat("Save noisy CN profiles (base=truth) in WIG format...\n")
        sample_cell_ID <- simulation$sample$sample_cell_ID
        noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
        for (cell in 1:length(sample_cell_ID)) {
            cell_ID <- sample_cell_ID[cell]
            filename <- paste(model_prefix_tmp, "/", model_prefix_tmp, "_noisy_cn_profiles_long_", iteration, "_", cell_ID, ".wig", sep = "")
            p0_write_cn_as_wig(filename, noisy_cn_profiles_long, cell_ID)
        }
    }
    #-----------Save noisy CN profiles on base of neuvar profiles in WIG
    if ((HMM == TRUE) & (neutral_variations == TRUE) & ((model_readcount_base == "neuvar") | (model_readcount_base == "all"))) {
        if (report_progress == TRUE) cat("Save noisy CN profiles (base=neuvar) in WIG format...\n")
        sample_cell_ID <- simulation$neutral_variations$sample$sample_cell_ID
        noisy_cn_profiles_long <- simulation$neutral_variations$sample$noisy_cn_profiles_long
        for (cell in 1:length(sample_cell_ID)) {
            cell_ID <- sample_cell_ID[cell]
            filename <- paste(model_prefix_tmp, "/", model_prefix_tmp, "_noisy_neuvar_cn_profiles_long_", iteration, "_", cell_ID, ".wig", sep = "")
            p0_write_cn_as_wig(filename, noisy_cn_profiles_long, cell_ID)
        }
    }
    #------------------------------------Run HMMcopy on noisy readcounts
    if (HMM == TRUE) {
        if (report_progress == TRUE) cat("Run HMMcopy on noisy CN profiles of all bases...\n")
        #   Run HMMcopy for each individual cell
        p0_write_gc_map_as_wig(filename_gc = paste(model_prefix_tmp, "/", model_prefix_tmp, "_gc_", iteration, ".wig", sep = ""), filename_map = paste(model_prefix_tmp, "/", model_prefix_tmp, "_map_", iteration, ".wig", sep = ""))
        if (HMM_containner == "singularity") {
            system(paste("sh hmmcopy_singularity.sh ", model_prefix_tmp, " ", iteration, sep = ""))
        } else if (HMM_containner == "docker") {
            system(paste("sh hmmcopy_docker.sh ", model_prefix_tmp, " ", iteration, sep = ""))
        }
        #   Update the simulation with HMMcopy inferences
        simulation <- p0_append_with_hmm(
            simulation = simulation,
            model_tmp = model_prefix_tmp,
            iteration = iteration,
            pseudo_corrected_readcount = pseudo_corrected_readcount,
            model_readcount_base = model_readcount_base
        )
    }
    # ======================================OUTPUT FILES FROM SIMULATION
    #----------------------------------------Save the simulation package
    if (save_simulation == TRUE) {
        if (report_progress == TRUE) cat("Save simulation package...\n")
        save(simulation, file = paste(folder_workplace_tmp, model_prefix_tmp, "_simulation_", iteration, ".rda", sep = ""))
    }
    #---------------------------Save the clonal identities of every cell
    if (save_cn_clones == TRUE) {
        if (report_progress == TRUE) cat("Save table of cell-clone mapping...\n")
        table_cell_clone <- simulation$sample$table_cell_clone
        filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_profiles_clonal_mapping_", iteration, ".csv", sep = "")
        write.csv(table_cell_clone, filename, row.names = FALSE)
    }
    #----------------------Save the table of CN events in cell phylogeny
    if (internal_nodes_cn_info == TRUE) {
        if (report_progress == TRUE) cat("Save table of CN events in cell phylogeny...\n")
        hclust_CN_events <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_CN_events
        filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_events_", iteration, ".csv", sep = "")
        write.csv(hclust_CN_events, filename, row.names = FALSE)
    }
    #-----------------------------Save the sampled cells' phylogeny tree
    if (save_newick_tree == TRUE) {
        if (report_progress == TRUE) cat("Save sampled cells' phylogeny in Newick format...\n")
        cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
        filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cell_phylogeny_", iteration, ".newick", sep = "")
        write(hc2Newick_MODIFIED(cell_phylogeny_hclust), file = filename)
        clone_phylogeny_hclust <- simulation$sample_phylogeny$clone_phylogeny_hclust
        filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_clone_phylogeny_", iteration, ".newick", sep = "")
        write(hc2Newick(clone_phylogeny_hclust), file = filename)
    }
    #---------------------------------------Save the sampled CN profiles
    if (save_cn_profile == TRUE) {
        if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save true CN profiles in long format...\n")
            cn_profiles_long <- simulation$sample$cn_profiles_long
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_profiles_long_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_long, filename, row.names = FALSE)
        }
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save true CN profiles in wide format...\n")
            cn_profiles_wide <- simulation$sample$cn_profiles_wide
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_profiles_wide_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_wide, filename, row.names = FALSE)
        }
    }
    #---------------Save the sampled CN profiles with neutral variations
    if ((save_cn_profile == TRUE) & (neutral_variations == TRUE)) {
        if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save true CN profiles in long format...\n")
            cn_profiles_long <- simulation$neutral_variations$sample$cn_profiles_long
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_profiles_long_neuvar_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_long, filename, row.names = FALSE)
        }
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save true CN profiles in wide format...\n")
            cn_profiles_wide <- simulation$neutral_variations$sample$cn_profiles_wide
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_cn_profiles_wide_neuvar_", iteration, ".csv", sep = "")
            write.csv(cn_profiles_wide, filename, row.names = FALSE)
        }
    }
    #--------------------Save noisy CN profiles on base of true profiles
    if ((save_cn_profile == TRUE) & (model_readcount == TRUE) & ((model_readcount_base == "truth") | (model_readcount_base == "all"))) {
        if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save noisy CN profiles (base=truth) in long format...\n")
            noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_noisy_cn_profiles_long_", iteration, ".csv", sep = "")
            write.csv(noisy_cn_profiles_long, filename, row.names = FALSE)
        }
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save noisy CN profiles (base=truth) in wide format...\n")
            noisy_cn_profiles_wide <- simulation$sample$noisy_cn_profiles_wide
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_noisy_cn_profiles_wide_", iteration, ".csv", sep = "")
            write.csv(noisy_cn_profiles_wide, filename, row.names = FALSE)
        }
    }
    #------------------Save noisy CN profiles on base of neuvar profiles
    if ((save_cn_profile == TRUE) & (neutral_variations == TRUE) & (model_readcount == TRUE) & ((model_readcount_base == "neuvar") | (model_readcount_base == "all"))) {
        if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save noisy CN profiles (base=neuvar) in long format...\n")
            noisy_cn_profiles_long <- simulation$neutral_variations$sample$noisy_cn_profiles_long
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_noisy_neuvar_cn_profiles_long_", iteration, ".csv", sep = "")
            write.csv(noisy_cn_profiles_long, filename, row.names = FALSE)
        }
        if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
            if (report_progress == TRUE) cat("Save noisy CN profiles (base=neuvar) in wide format...\n")
            noisy_cn_profiles_wide <- simulation$neutral_variations$sample$noisy_cn_profiles_wide
            filename <- paste(folder_workplace_tmp, model_prefix_tmp, "_noisy_neuvar_cn_profiles_wide_", iteration, ".csv", sep = "")
            write.csv(noisy_cn_profiles_wide, filename, row.names = FALSE)
        }
    }
    #---------------------------------------Return the simulation result
    print("HERE")
    if (length(output_variables) == 0) {
        return()
    } else if (any(output_variables == "all")) {
        return(simulation)
    } else {
        simulation_output <- vector("list", length(names(simulation)))
        names(simulation_output) <- names(simulation)
        for (i in 1:length(names(simulation))) {
            phase <- names(simulation)[i]
            vec_variables <- names(simulation[[phase]])
            locs <- which(vec_variables %in% output_variables)
            simulation_output[[phase]] <- simulation[[phase]][locs]
        }
        return(simulation_output)



        # simulation_output <- list()
        # if ("all_sample_genotype" %in% output_variables) simulation_output$sample$all_sample_genotype <- simulation$sample$all_sample_genotype
        # if ("sample_cell_ID" %in% output_variables) simulation_output$sample$sample_cell_ID <- simulation$sample$sample_cell_ID
        # if ("sample_genotype_unique" %in% output_variables) simulation_output$sample$sample_genotype_unique <- simulation$sample$sample_genotype_unique
        # if ("sample_genotype_unique_profile" %in% output_variables) simulation_output$sample$sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
        # if ("sample_genotype_unique_drivers" %in% output_variables) simulation_output$sample$sample_genotype_unique_drivers <- simulation$sample$sample_genotype_unique_drivers
        # if ("cell_phylogeny_hclust" %in% output_variables) simulation_output$sample_phylogeny$cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
        # if ("cn_profiles_long" %in% output_variables) {
        #     sample <- simulation$sample
        #     if (is.null(sample[["cn_profiles_long"]])) simulation <- p2_cn_profiles_long(simulation)
        #     simulation_output$sample$cn_profiles_long <- simulation$sample$cn_profiles_long
        # }
        # if ("cn_profiles_wide" %in% output_variables) {
        #     sample <- simulation$sample
        #     if (is.null(sample[["cn_profiles_wide"]])) simulation <- p2_cn_profiles_wide(simulation)
        #     simulation_output$sample$cn_profiles_wide <- simulation$sample$cn_profiles_wide
        # }
        # return(simulation_output)
    }
}
