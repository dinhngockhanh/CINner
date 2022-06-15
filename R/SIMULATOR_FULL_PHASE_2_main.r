# =============================================PHASE 2: SAMPLE PHYLOGENY
SIMULATOR_FULL_PHASE_2_main <- function(package_clonal_evolution, report_progress) {
    #-----------------------------------------Input the clonal evolution
    T_current <- package_clonal_evolution$T_current
    N_cells_current <- package_clonal_evolution$N_cells_current
    N_clones <- package_clonal_evolution$N_clones

    evolution_traj_time <- package_clonal_evolution$evolution_traj_time
    evolution_traj_clonal_ID <- package_clonal_evolution$evolution_traj_clonal_ID
    evolution_traj_population <- package_clonal_evolution$evolution_traj_population

    for (row in 1:nrow(Table_sampling)) {
        if (Table_sampling$Cell_count[row] == Inf) {
            loc <- which.min(abs(evolution_traj_time - Table_sampling$T_sample[row]))
            N_cells_total <- sum(evolution_traj_population[[loc]])
            Table_sampling$Cell_count[row] <- N_cells_total
        }
    }
    #---------------------------Find a random sample of final population
    Table_sampling$T_sample_real <- Table_sampling$T_sample

    all_sample_genotype <- c()
    all_sample_ID <- c()
    all_sample_sampled_time <- c()
    for (sample in 1:nrow(Table_sampling)) {
        N_sample <- Table_sampling$Cell_count[sample]
        ID_sample <- Table_sampling$Sample_ID[sample]
        T_sample <- Table_sampling$T_sample[sample]



        loc <- which.min(abs(evolution_traj_time - T_sample))
        Table_sampling$T_sample_real[sample] <- evolution_traj_time[loc]



        vec_clonal_ID <- evolution_traj_clonal_ID[[loc]]
        vec_clonal_population <- evolution_traj_population[[loc]]
        vec_population <- c()
        for (i in 1:length(vec_clonal_ID)) {
            clone <- vec_clonal_ID[i]
            clonal_population <- vec_clonal_population[i]
            vec_population <- c(vec_population, rep(clone, 1, clonal_population))
        }

        sample_genotype <- sample(x = vec_population, size = N_sample, replace = FALSE)
        all_sample_genotype <- c(all_sample_genotype, sample_genotype)
        all_sample_ID <- c(all_sample_ID, rep(ID_sample, N_sample))
        all_sample_sampled_time <- c(all_sample_sampled_time, rep(T_sample, N_sample))

        if (report_progress == TRUE) {
            cat(paste("Detected ", length(unique(sample_genotype)), " clones in sample ", ID_sample, "\n", sep = ""))
        }
    }
    if (report_progress == TRUE) {
        cat(paste("Detected ", length(unique(all_sample_genotype)), " clones in all samples\n", sep = ""))
    }
    #-----------------------------Create CN object for the sampled cells
    #---Find the CN profiles for each clone found in the sample
    sample_genotype_unique <- unique(all_sample_genotype)
    sample_genotype_unique_profile <- list()
    for (i_clone in 1:length(sample_genotype_unique)) {
        #   Extract CN information for the clone from clonal evolution data
        clone_ID <- sample_genotype_unique[i_clone]
        genotype_unique_profile <- get_cn_profile(package_clonal_evolution, clone_ID)
        sample_genotype_unique_profile[[i_clone]] <- genotype_unique_profile
    }
    #---Find the CN profiles for each cell in the sample
    sample_cell_ID <- c()
    sample_clone_ID <- all_sample_genotype
    sample_time <- all_sample_sampled_time
    #
    for (i_cell in 1:length(all_sample_genotype)) {
        sample_ID <- all_sample_ID[i_cell]
        cell_ID <- paste(sample_ID, "-Library-", as.character(i_cell), "-", as.character(i_cell), sep = "")
        sample_cell_ID[i_cell] <- cell_ID
    }
    #----------------------Give each clone a character index (A,B,C,...)
    sample_clone_ID_numeric <- all_sample_genotype
    sample_clone_ID_unique_numeric <- unique(sample_clone_ID_numeric)
    sample_clone_ID_unique_letters <- rep("", length = length(sample_clone_ID_unique_numeric))
    for (i in 1:length(sample_clone_ID_unique_letters)) {
        if (i <= 26) {
            sample_clone_ID_unique_letters[i] <- LETTERS[as.numeric(i)]
        } else {
            if ((i %% 26) == 0) {
                sample_clone_ID_unique_letters[i] <- paste(LETTERS[as.numeric(floor((i - 1) / 26))], "Z", sep = "")
            } else {
                sample_clone_ID_unique_letters[i] <- paste(LETTERS[as.numeric(floor((i - 1) / 26))], LETTERS[as.numeric(i %% 26)], sep = "")
            }
        }
    }


    table_clone_ID_vs_letters <- data.frame(sample_clone_ID_unique_numeric, sample_clone_ID_unique_letters)
    colnames(table_clone_ID_vs_letters) <- c("Clone_ID_number", "Clone_ID_letter")
    sample_clone_ID_letters <- c()
    for (i_clone in 1:length(sample_clone_ID_unique_numeric)) {
        clone_ID_numeric <- sample_clone_ID_unique_numeric[i_clone]
        clone_ID_letters <- sample_clone_ID_unique_letters[i_clone]
        vec_cell_ID <- which(sample_clone_ID_numeric == clone_ID_numeric)
        sample_clone_ID_letters[vec_cell_ID] <- clone_ID_letters
    }
    #-----------------------------Output package of data from simulation
    output <- list()
    output$sample_cell_ID <- sample_cell_ID
    output$sample_clone_ID <- sample_clone_ID
    output$sample_clone_ID_letters <- sample_clone_ID_letters
    output$table_clone_ID_vs_letters <- table_clone_ID_vs_letters
    output$sample_time <- sample_time

    output$all_sample_genotype <- all_sample_genotype
    output$all_sample_sampled_time <- all_sample_sampled_time
    output$all_sample_ID <- all_sample_ID
    output$sample_genotype_unique <- sample_genotype_unique
    output$sample_genotype_unique_profile <- sample_genotype_unique_profile

    output$Table_sampling <- Table_sampling

    return(output)
}
