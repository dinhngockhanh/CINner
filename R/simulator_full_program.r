simulator_full_program <- function(model = "",
                                   n_simulations = 0,
                                   stage_final = 0,
                                   n_clones_min = 0,
                                   n_clones_max = Inf,
                                   save_cn_profile = FALSE,
                                   format_cn_profile = "both",
                                   save_newick_tree = FALSE) {
    for (i in 1:n_simulations) {
        cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        cat(paste("BEGINNING SIMULATION-", i, "...\n", sep = ""))
        simulation <- one_simulation(model, stage_final, save_cn_profile, format_cn_profile)
        #------------------------------------Save the simulation package
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        save(simulation, file = filename)
        #-----------------------------------Save the sampled CN profiles
        if (save_cn_profile == TRUE) {
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                cn_profiles_long <-
                    simulation$sample$cn_profiles_long
                filename <- paste(model, "_cn_profiles_long_",
                    i, ".csv",
                    sep = ""
                )
                write.csv(cn_profiles_long,
                    filename,
                    row.names = FALSE
                )
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                cn_profiles_wide <-
                    simulation$sample$cn_profiles_wide
                filename <- paste(model, "_cn_profiles_wide_",
                    i, ".csv",
                    sep = ""
                )
                write.csv(cn_profiles_wide,
                    filename,
                    row.names = FALSE
                )
            }
        }
        #-------------------------Save the sampled cells' phylogeny tree
        if (save_newick_tree == TRUE) {
            cell_phylogeny_hclust <-
                simulation$sample_phylogeny$cell_phylogeny_hclust

            filename <- paste(model, "_cell_phylogeny_", i, ".newick", sep = "")
            write(hc2Newick(cell_phylogeny_hclust), file = filename)
        }
    }
}

one_simulation <- function(model, stage_final, save_cn_profile, format_cn_profile) {
    #-------------------------------------------Load model variables
    SIMULATOR_VARIABLES_for_simulation(model)
    #-----------------------------------------Produce one simulation
    flag_success <- 0
    while (flag_success == 0) {
        #------------------------------Simulate the clonal evolution
        cat("\nStage 1: clonal evolution...\n")
        output <- SIMULATOR_FULL_PHASE_1_main()
        flag_success <- output$flag_success
        package_clonal_evolution <- output$package_clonal_evolution
        if (flag_success == 0) {
            print("SIMULATION CONDITIONS NOT SATISFIED; REDOING...")
            next
        }
        #----------------------------------------Simulate the sample
        if (stage_final >= 2) {
            cat("\nStage 2: sampling...\n")
            package_sample <-
                SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
            N_clones <- nrow(package_sample$table_clone_ID_vs_letters)
            if ((N_clones < n_clones_min) | (N_clones > n_clones_max)) {
                flag_success <- 0
                cat("SIMULATION CONDITIONS NOT SATISFIED; REDOING...\n")
                next
            }
        }
        #-----------------------Simulate the phylogeny of the sample
        if (stage_final >= 3) {
            cat("\nStage 3: sample phylogeny...\n")
            package_sample_phylogeny <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution, package_sample)
        }
    }
    #--------------------------Save the simulation to output package
    if (stage_final >= 1) {
        simulation <- list()
        simulation$clonal_evolution <- package_clonal_evolution
    }
    if (stage_final >= 2) {
        simulation$sample <- package_sample
        if (save_cn_profile == TRUE) {
            if ((format_cn_profile == "long") | (format_cn_profile == "both")) {
                cat("\nExtra: build CN profile table in long format...\n")
                simulation <- p2_cn_profiles_long(simulation)
            }
            if ((format_cn_profile == "wide") | (format_cn_profile == "both")) {
                cat("\nExtra: build CN profile table in wide format...\n")
                simulation <- p2_cn_profiles_wide(simulation)
            }
        }
    }
    if (stage_final >= 3) {
        simulation$sample_phylogeny <- package_sample_phylogeny
    }
    return(simulation)
}
