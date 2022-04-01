# ==================================================CREATE ONE SIMULATION
SIMULATOR_FULL_PROGRAM_one_simulation <- function(model = "",
                                                  stage_final = 0,
                                                  N_clones_min = 0,
                                                  N_clones_max = Inf) {
    SIMULATOR_VARIABLES_for_simulation(model)

    flag_success <- 0
    while (flag_success == 0) {
        #------------------------------------------Simulate the clonal evolution
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("CLONAL EVOLUTION...")
        start.time <- Sys.time()
        output <- SIMULATOR_FULL_PHASE_1_main()
        flag_success <- output[[1]]
        package_clonal_evolution <- output[[2]]
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
        if (flag_success == 0) {
            print("SIMULATION CONDITION NOT SATISFIED; REDOING...")
            next
        }
        #----------------------------------------------------Simulate the sample
        if (stage_final >= 2) {
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("SAMPLING...")
            start.time <- Sys.time()
            package_sample <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
            N_clones <- nrow(package_sample[[5]])
            if ((N_clones < N_clones_min) | (N_clones > N_clones_max)) {
                flag_success <- 0
            }
        }
        if (flag_success == 0) {
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("SIMULATION CONDITION NOT SATISFIED; REDOING...")
            next
        }
        #-----------------------------------Simulate the phylogeny of the sample
        if (stage_final >= 3) {
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("SAMPLE PHYLOGENY...")
            start.time <- Sys.time()
            package_sample_phylogeny <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution, package_sample)
            end.time <- Sys.time()
            time.taken <- end.time - start.time
            print(time.taken)
        }
    }
    #----------------------------------Save the simulation to output package
    if (stage_final >= 1) {
        package_output <- list()
        package_output[[1]] <- package_clonal_evolution
    }
    if (stage_final >= 2) {
        package_output[[2]] <- package_sample
    }
    if (stage_final >= 3) {
        package_output[[3]] <- package_sample_phylogeny
    }
    #-------------------------------------------Save the simulation to files
    #   Save the sample's cell CN profiles
    if (stage_final >= 2) {
        # ??????????????????????????????????????????????????????????????
        # ??????????????????????????????????????????????????????????????
        # ??????????????????????????????????????????????????????????????
        all_sample_genotype <- package_sample[[7]]
        all_sample_sampled_time <- package_sample[[8]]
        all_sample_ID <- package_sample[[9]]
        sample_genotype_unique <- package_sample[[10]]
        sample_genotype_unique_profile <- package_sample[[11]]
        #---Find the CN profiles for each cell in the sample
        sample_cell_ID <- c()
        sample_clone_ID <- all_sample_genotype
        sample_time <- all_sample_sampled_time

        for (i_cell in 1:length(all_sample_genotype)) {
            if (i_cell %% 10 == 0) {
                print(paste("Cell-", i_cell, sep = ""))
            }

            sample_ID <- all_sample_ID[i_cell]
            clone_ID <- sample_clone_ID[i_cell]
            i_clone <- which(sample_genotype_unique == clone_ID)[1]
            #       Find the CN profile for this cell
            cell_genotype_profile <- sample_genotype_unique_profile[[i_clone]]
            #       Add column for cell ID
            cell_ID <- paste(sample_ID, "-Library-", as.character(i_cell), "-", as.character(i_cell), sep = "")
            sample_cell_ID[i_cell] <- cell_ID
            cell_genotype_profile <- cbind(cell_genotype_profile, rep(cell_ID, nrow(cell_genotype_profile)))
            names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
            #       Update table of CN profiles for all cells in the sample
            if (i_cell == 1) {
                sample_genotype_profile <- cell_genotype_profile
            } else {
                sample_genotype_profile <- rbind(sample_genotype_profile, cell_genotype_profile)
            }
        }

        package_sample[[1]] <- sample_genotype_profile
        package_output[[2]] <- package_sample
        # ??????????????????????????????????????????????????????????????
        # ??????????????????????????????????????????????????????????????
        # ??????????????????????????????????????????????????????????????


        # sample_genotype_profile <- package_sample[[1]]
        filename <- paste(model, "-output-cn-profiles", ".rda", sep = "")
        save(sample_genotype_profile, file = filename)
    }
    #   Save the sample's cell phylogeny
    if (stage_final >= 3) {
        phylogeny_clustering_truth <- package_sample_phylogeny[[1]]
        filename <- paste(model, "-output-clustering", ".rda", sep = "")
        save(phylogeny_clustering_truth, file = filename)
    }
    #------------------------------------------Output the simulation package
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("DONE WITH SIMULATION...")
    return(package_output)
}
