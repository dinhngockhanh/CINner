#' @export
fitting_bulk_gene <- function(model_name,
                              model_variables,
                              driver_DATA,
                              list_parameters,
                              list_targets,
                              ABC_simcount = 10000,
                              n_cores = NULL,
                              n_samples = NULL,
                              R_libPaths = NULL) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(dplyr)
    if (is.null(n_cores)) {
        n_cores <- detectCores() - 1
    } else {
        n_cores <- n_cores - 1
    }
    numCores <- n_cores + 1
    if (is.null(n_samples) == TRUE) {
        n_samples <- 100
    }
    #---------------------------Define functions for the fitting routine
    #---Function to assign parameters to proper positions
    assign_paras <- function(model_variables, parameter_IDs, parameters) {
        for (i in 1:length(parameter_IDs)) {
            parameter_ID <- parameter_IDs[i]
            if (parameter_ID %in% model_variables$general_variables$Variable) {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameters[i]
            } else if (parameter_ID %in% model_variables$driver_library$Gene_ID) {
                model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == parameter_ID)] <- parameters[i]
            }
            model_variables <- BUILD_driver_library(
                model_variables = model_variables,
                table_gene_selection_rates = model_variables$driver_library
            )
        }
        return(model_variables)
    }
    #---Function to extract mut/amp/del freqs from each simulation
    get_gene_freqs <- function(SIMS_raw, list_targets, driver_library) {
        count_mut <- matrix(0, nrow = length(SIMS_raw), ncol = length(list_targets))
        count_gain <- matrix(0, nrow = length(SIMS_raw), ncol = length(list_targets))
        count_loss <- matrix(0, nrow = length(SIMS_raw), ncol = length(list_targets))
        #------Find mut/gain/loss status of each gene in each simulation
        for (iteration in 1:length(SIMS_raw)) {
            simulation <- SIMS_raw[[iteration]]
            all_sample_genotype <- simulation$sample$all_sample_genotype
            sample_genotype_unique <- simulation$sample$sample_genotype_unique
            sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
            sample_genotype_unique_drivers <- simulation$sample$sample_genotype_unique_drivers
            #   Get the genotype with the highest clonal percentage in sample
            tmp <- as.data.frame(table(all_sample_genotype))
            max_freq <- max(tmp$Freq)
            vec_loc <- which(tmp$Freq == max_freq)
            if (length(vec_loc) > 1) {
                loc <- sample(vec_loc, 1)
            } else {
                loc <- vec_loc
            }
            max_genotype <- tmp$all_sample_genotype[loc]
            #   Get this clone's mutation profile
            max_genotype_drivers <- sample_genotype_unique_drivers[[max_genotype]]
            #   Get this clone's CN profile
            max_genotype_CN <- sample_genotype_unique_profile[[max_genotype]]
            #   Get this clone's ploidy
            max_genotype_ploidy <- max(1, round(mean(max_genotype_CN$copy)))
            #   Update counts of gains/losses
            for (i in 1:length(list_targets)) {
                Gene_ID <- list_targets[i]
                chr <- driver_library$Chromosome[which(driver_library$Gene_ID == Gene_ID)]
                bin <- driver_library$Bin[which(driver_library$Gene_ID == Gene_ID)]
                Gene_CN <- max_genotype_CN$copy[which(max_genotype_CN$chr == chr)][bin]
                if (Gene_CN > max_genotype_ploidy) {
                    count_gain[iteration, i] <- 1
                } else if (Gene_CN < max_genotype_ploidy) {
                    count_loss[iteration, i] <- 1
                }
            }
            #   Update counts of mutations
            if (nrow(max_genotype_drivers) >= 1) {
                for (i in 1:nrow(max_genotype_drivers)) {
                    tmp <- max_genotype_drivers[i, 1]
                    Gene_ID <- driver_library$Gene_ID[tmp]
                    count_mut[iteration, which(list_targets == Gene_ID)] <- 1
                }
            }
        }
        #----Output mut/gain/loss status of each gene in each simulation
        sim_freqs <- list()
        sim_freqs$count_mut <- count_mut
        sim_freqs$count_gain <- count_gain
        sim_freqs$count_loss <- count_loss
        return(sim_freqs)
    }
    #---Objective function for ABC fitting
    func_ABC <- function(parameters, parameter_IDs, model_variables, list_targets) {
        #   Assign parameters in model variables
        model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
        #   Get driver gene library
        driver_library <- model_variables$driver_library
        #   Make simulations
        SIMS_raw <- simulator_full_program(
            model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
            save_simulation = FALSE, report_progress = TRUE,
            output_variables = c("all_sample_genotype", "sample_genotype_unique", "sample_genotype_unique_profile", "sample_genotype_unique_drivers"),
            build_cn = FALSE
        )
        #   Statistics = frequencies of mut/amp/del of each driver gene
        sim_freqs <- get_gene_freqs(SIMS_raw, list_targets, driver_library)
        count_mut <- sim_freqs$count_mut
        count_gain <- sim_freqs$count_gain
        count_loss <- sim_freqs$count_loss
        freq_mut <- 100 * colMeans(count_mut)
        freq_gain <- 100 * colMeans(count_gain)
        freq_loss <- 100 * colMeans(count_loss)
        stat <- c(freq_mut, freq_gain, freq_loss)
    }
    #---Function to choose one best parameter from a posterior
    get_best_para <- function(data_rf, model_rf, obs_rf, post_rf) {
        df_dist <- densityPlot_df(model_rf, obs_rf, data_rf)
        best_para <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
    }
    #---------------------------------List of parameter IDs to be fitted
    parameter_IDs <- list_parameters$Variable
    #------------------------Find target statistics for all driver genes
    #   Get the target statistics (= mut/gain/loss freqs of each gene)
    DATA_target <- c(driver_DATA$Freq_mut, driver_DATA$Freq_gain, driver_DATA$Freq_loss)



    # parameters <- c(1e-4, 1e-4, 1.5, rep(1.1, length(list_targets) - 1))
    # tmp <- func_ABC(parameters, parameter_IDs, model_variables, list_targets)
    # print(matrix(tmp, nrow = 3, byrow = TRUE))
    # print(matrix(DATA_target, nrow = 3, byrow = TRUE))
    # return()



    # # =============================================CREATE REFERENCE DATA
    # #---------------------------------------Simulate table of parameters
    # sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    # for (col in 1:ncol(sim_param)) {
    #     sim_param[, col] <- runif(ABC_simcount, min = as.numeric(list_parameters$Lower_bound[col]), max = as.numeric(list_parameters$Upper_bound[col]))
    # }
    # #-----------------------------------------------Make reference table
    # start_time <- Sys.time()
    # #   Configure parallel pool
    # cl <- makePSOCKcluster(numCores - 1)
    # cat(paste("\nParallel cluster with ", numCores - 1, " cores...\n", sep = ""))
    # cat("Creating reference table for ABC...\n")
    # if (is.null(R_libPaths) == FALSE) {
    #     R_libPaths <<- R_libPaths
    #     clusterExport(cl, varlist = c("R_libPaths"))
    #     clusterEvalQ(cl = cl, .libPaths(R_libPaths))
    # }
    # #   Prepare input parameters for creating reference data
    # sim_param <<- sim_param
    # parameter_IDs <<- parameter_IDs
    # model_variables <<- model_variables
    # list_targets <<- list_targets
    # func_ABC <<- func_ABC
    # assign_paras <<- assign_paras
    # clusterExport(cl, varlist = c(
    #     "sim_param", "parameter_IDs", "model_variables", "list_targets", "func_ABC", "assign_paras",
    #     "BUILD_driver_library", "simulator_full_program", "one_simulation",
    #     "SIMULATOR_VARIABLES_for_simulation",
    #     "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
    #     "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
    #     "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
    #     "SIMULATOR_FULL_PHASE_2_main",
    #     "SIMULATOR_FULL_PHASE_3_main",
    #     "get_cn_profile", "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model"
    # ))
    # e <- new.env()
    # e$libs <- .libPaths()
    # clusterExport(cl, "libs", envir = e)
    # clusterEvalQ(cl, .libPaths(libs))
    # #   Create simulated statistics in parallel
    # pbo <- pboptions(type = "txt")
    # sim_results_list <- pblapply(cl = cl, X = 1:ABC_simcount, FUN = function(iteration) {
    #     parameters <- sim_param[iteration, ]
    #     stat <- func_ABC(parameters, parameter_IDs, model_variables, list_targets)
    #     return(stat)
    # })
    # stopCluster(cl)
    # #   Group simulated statistics into one table
    # sim_stat <- matrix(0, nrow = ABC_simcount, ncol = length(DATA_target))
    # for (row in 1:ABC_simcount) {
    #     stat <- sim_results_list[[row]]
    #     sim_stat[row, ] <- stat
    # }
    # end_time <- Sys.time()
    # print(end_time - start_time)
    # #---------------------------Save the parameters and their statistics
    # ABC_input <- list()
    # ABC_input$model_variables <- model_variables
    # ABC_input$n_samples <- n_samples
    # ABC_input$parameter_IDs <- parameter_IDs
    # ABC_input$sim_param <- sim_param
    # ABC_input$sim_stat <- sim_stat
    # ABC_input$DATA_target <- DATA_target
    # filename <- paste(model_name, "_ABC_input.rda", sep = "")
    # save(ABC_input, file = filename)



    # # ====================================FITTING WITH ABC RANDOM FOREST
    # filename <- paste(model_name, "_ABC_input.rda", sep = "")
    # load(filename)
    # model_variables <- ABC_input$model_variables
    # n_samples <- ABC_input$n_samples
    # parameter_IDs <- ABC_input$parameter_IDs
    # sim_param <- ABC_input$sim_param
    # sim_stat <- ABC_input$sim_stat
    # #--------------------------------------------Fit parameters with ABC
    # #   Dataframe for parameters for reference
    # all_paras <- data.frame(sim_param)
    # colnames(all_paras) <- parameter_IDs
    # #   Dataframe for corresponding statistics for reference
    # all_data <- data.frame(sim_stat)
    # colnames(all_data) <- c(paste0("mut_", list_targets), paste0("gain_", list_targets), paste0("loss_", list_targets))
    # #   Dataframe for data observation
    # obs_rf <- data.frame(matrix(DATA_target, nrow = 1))
    # colnames(obs_rf) <- c(paste0("mut_", list_targets), paste0("gain_", list_targets), paste0("loss_", list_targets))
    # #   Fit each parameter with ABC-rf
    # layout <- matrix(NA, nrow = 7, ncol = ceiling(length(parameter_IDs) / 7))
    # gs <- list()
    # id <- 0
    # for (para in 1:length(parameter_IDs)) {
    #     para_ID <- parameter_IDs[para]
    #     cat(paste("ABC for parameter ", para_ID, " [", para, "/", length(parameter_IDs), "]", "\n", sep = ""))
    #     #   Train the random forest
    #     data_rf <- cbind(all_paras[para_ID], all_data)
    #     colnames(data_rf)[1] <- "para"
    #     f <- as.formula("para ~.")
    #     model_rf <- regAbcrf(formula = f, data_rf, paral = TRUE, ncores = n_cores)
    #     #   Predict posterior distribution based on found random forest
    #     post_rf <- predict(model_rf, obs_rf, data_rf, paral = TRUE, ncores = n_cores)
    #     #   Choose best value from posterior distribution
    #     best_rf <- get_best_para(data_rf, model_rf, obs_rf, post_rf)
    #     #   Save results for fitting this parameter
    #     ABC_output <- list()
    #     ABC_output$para_ID <- para_ID
    #     ABC_output$data_rf <- data_rf
    #     ABC_output$model_rf <- model_rf
    #     ABC_output$obs_rf <- obs_rf
    #     ABC_output$post_rf <- post_rf
    #     ABC_output$best_rf <- best_rf
    #     filename <- paste(model_name, "_ABC_output_", para_ID, ".rda", sep = "")
    #     save(ABC_output, file = filename)
    #     #   Plot the prior, posterior and chosen best parameter
    #     filename <- paste0(model_name, "_ABC_", para_ID, ".jpeg")
    #     jpeg(filename, width = 2000, height = 1100)
    #     p <- densityPlot_MODIFIED(
    #         model_rf, obs_rf, data_rf,
    #         protocol = "arm",
    #         chosen_para = best_rf,
    #         color_prior = "lightblue", color_posterior = "darkblue", color_vline = "blue",
    #         main = para_ID
    #     )
    #     print(p)
    #     dev.off()
    #     #   Plot the prior, posterior and chosen best parameter for all variables
    #     id <- id + 1
    #     row <- id %% 7
    #     if (row == 0) row <- 7
    #     col <- ceiling(id / 7)
    #     layout[row, col] <- id
    #     gs[[id]] <- densityPlot_MODIFIED(
    #         model_rf, obs_rf, data_rf,
    #         protocol = "arm",
    #         fontsize = 20,
    #         chosen_para = best_rf,
    #         color_prior = "lightblue", color_posterior = "darkblue", color_vline = "blue",
    #         main = para_ID
    #     )
    # }
    # #   Plot the prior, posterior and chosen best parameter for all variables
    # filename <- paste0(model_name, "_ABC_all.jpeg")
    # jpeg(filename, width = 3000, height = 1500)
    # p <- grid.arrange(grobs = gs, layout_matrix = layout)
    # print(p)
    # dev.off()
    # =======================================ANALYSIS OF FITTING RESULTS
    #------------------Choose the best parameter set from all posteriors
    parameters_best <- rep(0, length(parameter_IDs))
    for (para in 1:length(parameter_IDs)) {
        para_ID <- parameter_IDs[para]
        filename <- paste(model_name, "_ABC_output_", para_ID, ".rda", sep = "")
        load(filename)
        best_rf <- ABC_output$best_rf
        parameters_best[para] <- best_rf
    }
    #-----------------Analysis of fitted driver frequencies against data
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters_best)
    #   Make simulations using best parameters
    SIMS_raw <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = 1000, stage_final = 2,
        # model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE,
        output_variables = c("all_sample_genotype", "sample_genotype_unique", "sample_genotype_unique_profile", "sample_genotype_unique_drivers"),
        build_cn = FALSE,
        compute_parallel = TRUE,
        R_libPaths = R_libPaths
    )
    sim_freqs <- get_gene_freqs(SIMS_raw, list_targets, model_variables$driver_library)
    count_mut <- sim_freqs$count_mut
    count_gain <- sim_freqs$count_gain
    count_loss <- sim_freqs$count_loss
    #   Compute averages of frequencies from simulations
    mean_freqs <- c(100 * colMeans(count_mut), 100 * colMeans(count_gain), 100 * colMeans(count_loss))
    #   Compute standard errors of frequencies from simulations
    N_bootstrap <- 10000
    stat_bootstrap <- matrix(0, nrow = N_bootstrap, ncol = length(mean_freqs))
    for (i in 1:N_bootstrap) {
        locs <- sample(1:length(SIMS_raw), length(SIMS_raw), replace = TRUE)
        count_mut_tmp <- count_mut[locs, ]
        count_gain_tmp <- count_gain[locs, ]
        count_loss_tmp <- count_loss[locs, ]
        stat_bootstrap[i, ] <- c(100 * colMeans(count_mut_tmp), 100 * colMeans(count_gain_tmp), 100 * colMeans(count_loss_tmp))
    }
    sd_freqs <- apply(stat_bootstrap, 2, sd)
    #   Make dataframe of data frequencies and simulated frequencies
    Gene_ID <- driver_DATA$Gene_ID
    df_freq <- data.frame(
        Gene_ID = rep(driver_DATA$Gene_ID, 6),
        Freq = c(DATA_target, mean_freqs),
        Type = rep(c(rep("Mutation", length(driver_DATA$Gene_ID)), rep("Gain", length(driver_DATA$Gene_ID)), rep("Loss", length(driver_DATA$Gene_ID))), 2),
        Source = c(rep("Data", length(DATA_target)), rep("Simulation", length(mean_freqs))),
        SD = c(rep(0, length(DATA_target)), sd_freqs)
    )
    #   Plot comparison vs data
    filename <- paste0(model_name, "_against_data.jpeg")
    jpeg(file = filename, width = 2000, height = 1100)
    p <- ggplot(df_freq, aes(x = Gene_ID, y = Freq, fill = Type, color = Type)) +
        geom_bar(stat = "identity", width = 0.5, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Simulation"), alpha = 0.6) +
        geom_errorbar(aes(ymin = Freq - SD, ymax = Freq + SD), width = 0.5, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Simulation")) +
        geom_point(size = 5, shape = 21, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Data")) +
        xlab("Gene") +
        ylab("Frequency") +
        scale_fill_manual(
            values = c(
                "Mutation" = "olivedrab",
                "Gain" = "#ff0000",
                "Loss" = "#0072B2"
            )
        ) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()

    #   ............................................................
    #   ............................................................
    #   ............................................................
    #   ............................................................
    #   ............................................................
    #   ............................................................
    #   ............................................................
}
