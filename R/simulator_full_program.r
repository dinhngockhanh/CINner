simulator_full_program <- function(model = "",
                                   n_simulations = 0,
                                   stage_final = 0,
                                   n_clones_min = 0,
                                   n_clones_max = Inf,
                                   save_cn_profile = FALSE,
                                   save_newick_tree = FALSE) {
    for (i in 1:n_simulations) {
        print("=======================================================")
        print(paste("BEGINNING SIMULATION-", i, "...", sep = ""))
        simulation <- one_simulation(model, stage_final)
        #------------------------------------Save the simulation package
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        save(simulation, file = filename)
        #-----------------------------------Save the sampled CN profiles
        if (save_cn_profile == TRUE) {
            sample_genotype_profiles <- simulation$sample$sample_genotype_profiles

            filename <- paste(model, "_cn_profiles_", i, ".csv", sep = "")
            write.csv(sample_genotype_profiles, filename, row.names = FALSE)

            # filename <- paste(model, "_cn_profiles_", i, ".rda", sep = "")
            # save(sample_genotype_profiles, file = filename)
        }
        #-------------------------Save the sampled cells' phylogeny tree
        if (save_newick_tree == TRUE) {
            cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust

            filename <- paste(model, "_cell_phylogeny_", i, ".newick", sep = "")
            write(hc2Newick(cell_phylogeny_hclust), file = filename)
        }
    }
}

one_simulation <- function(model, stage_final) {
    #-------------------------------------------Load model variables
    SIMULATOR_VARIABLES_for_simulation(model)
    #-----------------------------------------Produce one simulation
    flag_success <- 0
    while (flag_success == 0) {
        #------------------------------Simulate the clonal evolution
        print("Stage 1: clonal evolution...")
        output <- SIMULATOR_FULL_PHASE_1_main()
        flag_success <- output$flag_success
        package_clonal_evolution <- output$package_clonal_evolution
        if (flag_success == 0) {
            print("SIMULATION CONDITIONS NOT SATISFIED; REDOING...")
            next
        }
        #----------------------------------------Simulate the sample
        if (stage_final >= 2) {
            print("Stage 2: sampling...")
            package_sample <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
            N_clones <- nrow(package_sample$table_clone_ID_vs_letters)
            if ((N_clones < n_clones_min) | (N_clones > n_clones_max)) {
                flag_success <- 0
                print("SIMULATION CONDITIONS NOT SATISFIED; REDOING...")
                next
            }
        }
        #-----------------------Simulate the phylogeny of the sample
        if (stage_final >= 3) {
            print("Stage 3: sample phylogeny...")
            package_sample_phylogeny <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution, package_sample)
        }
    }
    #--------------------------Save the simulation to output package
    if (stage_final >= 1) {
        simulation <- list()
        simulation$clonal_evolution <- package_clonal_evolution
    }
    if (stage_final >= 2) {
        package_sample <- SIMULATOR_FULL_PHASE_2_copy_number_table(package_sample)
        simulation$sample <- package_sample
    }
    if (stage_final >= 3) {
        simulation$sample_phylogeny <- package_sample_phylogeny
    }
    return(simulation)
}
