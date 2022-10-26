#' @export
fitting_arm_PCAWG <- function(model_name,
                              model_variables,
                              copynumber_PCAWG,
                              list_parameters,
                              list_targets,
                              ABC_method = "rf",
                              ABC_simcount = 10000,
                              n_cores = NULL,
                              n_samples = NULL) {
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    if (is.null(n_samples) == TRUE) {
        n_samples <- length(unique(copynumber_PCAWG$donor_unique_id))
    }
    #------------------------------------------Get list of parameter IDs
    parameter_IDs <- list_parameters$Variable
    #------------------Function to assign parameters to proper positions
    assign_paras_PCAWG <- function(model_variables, parameter_IDs, parameters) {
        for (i in 1:length(parameter_IDs)) {
            parameter_ID <- parameter_IDs[i]
            if (parameter_ID %in% model_variables$general_variables$Variable) {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameters[i]
            } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)] <- parameters[i]
            }
        }
        return(model_variables)
    }
    #----------------------Function to extract arm-level gain/loss delta
    get_arm_gainloss <- function(bin_gainloss, copynumber_coordinates, list_targets) {
        arm_gain <- rep(0, length(list_targets))
        arm_loss <- rep(0, length(list_targets))
        for (i in 1:length(list_targets)) {
            chrom <- substr(list_targets[i], 1, nchar(list_targets[i]) - 1)
            arm <- substr(list_targets[i], nchar(list_targets[i]), nchar(list_targets[i]))
            vec_loc <- which(copynumber_coordinates$chr == chrom)
            delta_gain <- bin_gainloss$delta_gain[vec_loc]
            delta_loss <- bin_gainloss$delta_loss[vec_loc]
            if (arm == "p") {
                arm_gain[i] <- delta_gain[1]
                arm_loss[i] <- delta_loss[1]
            } else if (arm == "q") {
                arm_gain[i] <- delta_gain[length(delta_gain)]
                arm_loss[i] <- delta_loss[length(delta_loss)]
            }
        }
        arm_gainloss <- list()
        arm_gainloss$delta_gain <- arm_gain
        arm_gainloss$delta_loss <- arm_loss
        return(arm_gainloss)
    }
    #---------------------------------Objective function for ABC fitting
    func_ABC <- function(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets) {
        #   Assign parameters in model variables
        model_variables <- assign_paras_PCAWG(model_variables, parameter_IDs, parameters)
        #   Make simulations
        SIMS_chromosome <- simulator_full_program(
            model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
            save_simulation = FALSE, report_progress = TRUE,
            output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile")
        )
        #   Statistics = arm-level gain/loss delta
        SIMS_delta_genome_arms <- gainloss_SIMS(SIMS_chromosome, ploidy_normalization = FALSE)
        SIMS_arm_gainloss <- get_arm_gainloss(SIMS_delta_genome_arms, copynumber_coordinates, list_targets)
        stat <- c(SIMS_arm_gainloss$delta_gain, SIMS_arm_gainloss$delta_loss)
    }
    #-------------Function to choose one best parameter from a posterior
    get_best_para <- function(data_rf, model_rf, obs_rf, post_rf) {
        df_dist <- densityPlot_df(model_rf, obs_rf, data_rf)
        best_para <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
    }
    #---------------Find arm-level PCAWG gain/loss map for entire genome
    cn_info <- model_variables$cn_info
    #   Find genome coordinate from one simulation
    copynumber_sims <- simulator_full_program(
        model = model_variables,
        model_prefix = "TEST",
        n_simulations = 1,
        stage_final = 2,
        save_simulation = FALSE,
        report_progress = FALSE,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile")
    )
    simulation <- copynumber_sims[[1]]
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    CNbins_iteration <- sample_genotype_unique_profile[[1]]
    copynumber_coordinates <- CNbins_iteration[, 1:3]
    copynumber_coordinates$width <- copynumber_coordinates$end - copynumber_coordinates$start + 1
    #   Get gain/loss map from PCAWG for this chromosome at arm level
    PCAWG_delta_genome_arms <- gainloss_PCAWG(
        copynumber_PCAWG, copynumber_coordinates,
        ploidy_normalization = TRUE,
        use_rbindlist = TRUE,
        arm_level = TRUE, pos_centromeres = cn_info
    )
    #   Get the target statistics (= gain/loss on each chromosome arm)
    PCAWG_arm_gainloss <- get_arm_gainloss(PCAWG_delta_genome_arms, copynumber_coordinates, list_targets)
    PCAWG_target <- c(PCAWG_arm_gainloss$delta_gain, PCAWG_arm_gainloss$delta_loss)
    # ==============================================CREATE REFERENCE DATA
    # #---------------------------------------Simulate table of parameters
    # sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    # for (col in 1:ncol(sim_param)) {
    #     sim_param[, col] <- runif(ABC_simcount, min = as.numeric(list_parameters$Lower_bound[col]), max = as.numeric(list_parameters$Upper_bound[col]))
    # }
    # #-----------------------------------------------Make reference table
    # start_time <- Sys.time()
    # #   Configure parallel pool
    # cl <- makePSOCKcluster(n_cores)
    # cat("Creating reference table for ABC...\n")
    # sim_param <<- sim_param
    # parameter_IDs <<- parameter_IDs
    # model_variables <<- model_variables
    # gainloss_SIMS <<- gainloss_SIMS
    # func_ABC <<- func_ABC
    # assign_paras_PCAWG <<- assign_paras_PCAWG
    # clusterExport(cl, varlist = c(
    #     "sim_param", "parameter_IDs", "model_variables", "gainloss_SIMS", "func_ABC", "assign_paras_PCAWG",
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
    #     stat <- func_ABC(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets)
    #     return(stat)
    # })
    # stopCluster(cl)
    # #   Group simulated statistics into one table
    # sim_stat <- matrix(0, nrow = ABC_simcount, ncol = length(PCAWG_target))
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
    # ABC_input$PCAWG_target <- PCAWG_target
    # filename <- paste(model_name, "_ABC_input.rda", sep = "")
    # save(ABC_input, file = filename)
    # =====================================FITTING WITH ABC RANDOM FOREST
    # #--------------------------------------------Fit parameters with ABC
    # #   Dataframe for parameters for reference
    # all_paras <- data.frame(sim_param)
    # colnames(all_paras) <- parameter_IDs
    # #   Dataframe for corresponding statistics for reference
    # all_data <- data.frame(sim_stat)
    # colnames(all_data) <- paste("gainloss_", 1:ncol(all_data), sep = "")
    # #   Dataframe for PCAWG observation
    # obs_rf <- data.frame(matrix(PCAWG_target, nrow = 1))
    # colnames(obs_rf) <- paste("gainloss_", 1:ncol(obs_rf), sep = "")
    # #   Fit each parameter with ABC-rf
    # for (para in 1:length(parameter_IDs)) {
    #     para_ID <- parameter_IDs[para]
    #     cat(paste("ABC for parameter ", para_ID, "\n", sep = ""))
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
    #     filename <- paste("ABC_", para_ID, ".jpeg", sep = "")
    #     jpeg(filename, width = 2000, height = 1000)
    #     p <- densityPlot_MODIFIED(
    #         model_rf, obs_rf, data_rf,
    #         protocol = "arm",
    #         chosen_para = best_rf,
    #         color_prior = "lightblue", color_posterior = "darkblue", color_vline = "blue",
    #         main = para_ID
    #     )
    #     print(p)
    #     dev.off()
    # }
    # ========================================ANALYSIS OF FITTING RESULTS
    #------------------Choose the best parameter set from all posteriors
    parameters_best <- rep(0, length(parameter_IDs))
    for (para in 1:length(parameter_IDs)) {
        para_ID <- parameter_IDs[para]
        filename <- paste(model_name, "_ABC_output_", para_ID, ".rda", sep = "")
        load(filename)
        best_rf <- ABC_output$best_rf
        parameters_best[para] <- best_rf
    }
    # #-----------------------Analysis of fitted CN profiles against PCAWG
    # #   Assign parameters in model variables
    # model_variables <- assign_paras_PCAWG(model_variables, parameter_IDs, parameters_best)
    # #   Make simulations using best parameters
    # SIMS_chromosome <- simulator_full_program(
    #     model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
    #     save_simulation = FALSE, report_progress = TRUE, compute_parallel = TRUE,
    #     output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile")
    # )
    # #   Plot comparison vs bin-level PCAWG
    # filename <- paste(model_name, "_against_bin_PCAWG.jpeg", sep = "")
    # plot_gainloss(SIMS_chromosome, copynumber_PCAWG, filename, arm_level = FALSE, pos_centromeres = cn_info)
    # #   Plot comparison vs arm-level PCAWG
    # filename <- paste(model_name, "_against_arm_PCAWG.jpeg", sep = "")
    # plot_gainloss(SIMS_chromosome, copynumber_PCAWG, filename, arm_level = TRUE, pos_centromeres = cn_info)
    #--------------------------------------Analysis of fitted parameters
    #   Prepare dataframe of fitted selection rates and amp/del freqs
    plot_table <- data.frame(Arm = list_targets)
    plot_table$Amp_freq_spec <- 0
    plot_table$Del_freq_spec <- 0
    plot_table$Selection_rate <- 0
    for (row in 1:nrow(plot_table)) {
        arm <- plot_table$Arm[row]
        plot_table$Amp_freq_spec[row] <- PCAWG_target[which(list_targets == arm)]
        plot_table$Del_freq_spec[row] <- -PCAWG_target[which(list_targets == arm) + length(list_targets)]
        plot_table$Selection_rate[row] <- parameters_best[which(list_parameters == arm)]
    }

    print(plot_table)

    #   Configuration for subplots
    layout <- matrix(NA, nrow = 1, ncol = 2)
    gs <- list()
    plots <- list()
    id <- 0
    #------------Plot selection rates vs cancer-specific amp frequencies
    id <- id + 1
    layout[1, 1] <- id
    tmp <- cor.test(plot_table$Selection_rate, plot_table$Amp_freq_spec)
    Pearson_r <- tmp$estimate
    Pearson_p_val <- tmp$p.value
    grob1 <- grobTree(
        textGrob(
            paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
            x = 0.1, y = 0.9, just = "left",
            gp = gpar(col = "blue", fontsize = 30, fontface = "bold"),
        )
    )
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Amp_freq_spec)) +
        geom_point(size = 25, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        annotation_custom(grob1) +
        geom_text(aes(label = Arm), size = 10, color = "white") +
        xlab("Selection rates") +
        ylab(paste("Amplification frequencies (", model_name, ")", sep = "")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 30))
    gs[[id]] <- plots[[id]]
    #------------Plot selection rates vs cancer-specific del frequencies
    id <- id + 1
    layout[1, 2] <- id
    tmp <- cor.test(plot_table$Selection_rate, plot_table$Del_freq_spec)
    Pearson_r <- tmp$estimate
    Pearson_p_val <- tmp$p.value
    grob1 <- grobTree(
        textGrob(
            paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
            x = 0.9, y = 0.9, just = "right",
            gp = gpar(col = "blue", fontsize = 30, fontface = "bold")
        )
    )
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Del_freq_spec)) +
        geom_point(size = 25, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        annotation_custom(grob1) +
        geom_text(aes(label = Arm), size = 10, color = "white") +
        xlab("Selection rates") +
        ylab(paste("Deletion frequencies (", model_name, ")", sep = "")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 30))
    gs[[id]] <- plots[[id]]
    #----------------------------------------------------Print all plots
    filename <- paste(model_name, "_selection_rate_analysis.jpeg", sep = "")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()



    # #   Prepare ID lists to look up gain/loss freqs and fitted sel rates
    # list_targets_arm <- list_targets
    # list_targets_chrom <- list_targets
    # for (i in 1:length(list_targets_chrom)) {
    #     list_targets_chrom[i] <- substr(list_targets_chrom[i], 1, nchar(list_targets_chrom[i]) - 1)
    # }
    # list_parameters_arm <- parameter_IDs
    # list_parameters_chrom <- parameter_IDs
    # for (i in 1:length(list_parameters_chrom)) {
    #     list_parameters_chrom[i] <- substr(list_parameters_chrom[i], 1, nchar(list_parameters_chrom[i]) - 1)
    # }
    # #   Load the Davoli-Charm score table
    # CHARM_table <- data.frame(read_excel("Davoli_Charm_score.xlsx"))
    # CHARM_table <- CHARM_table[, c(1, 2, 3, 14, 15)]
    # colnames(CHARM_table) <- c("Arm", "Del_freq_all", "Amp_freq_all", "Charm.TSG.OG.score", "Charm.TSG.OG.Ess.score")
    # #   Load the Davoli-Chrom score table
    # CHROM_table <- data.frame(read_excel("Davoli_Chrom_score.xlsx"))
    # CHROM_table <- CHROM_table[, c(1, 2, 3, 14, 15)]
    # colnames(CHROM_table) <- c("Chromosome", "Del_freq_all", "Amp_freq_all", "Chrom.TSG.OG.score", "Chrom.TSG.OG.Ess.score")
    # #   Add cancer-specific del/amp frequencies to Davoli-Charm table
    # CHARM_table$Amp_freq_spec <- 0
    # CHARM_table$Del_freq_spec <- 0
    # for (row in 1:nrow(CHARM_table)) {
    #     arm <- CHARM_table$Arm[row]
    #     loc <- which(list_targets_arm == arm)
    #     CHARM_table$Amp_freq_spec[row] <- PCAWG_target[loc]
    #     CHARM_table$Del_freq_spec[row] <- -PCAWG_target[loc + length(list_targets_arm)]
    # }
    # #   Add cancer-specific del/amp frequencies to Davoli-Chrom table
    # CHROM_table$Amp_freq_spec <- 0
    # CHROM_table$Del_freq_spec <- 0
    # for (row in 1:nrow(CHROM_table)) {
    #     chrom <- CHROM_table$Chromosome[row]
    #     loc <- which(list_targets_chrom == chrom)
    #     CHROM_table$Amp_freq_spec[row] <- mean(PCAWG_target[loc])
    #     CHROM_table$Del_freq_spec[row] <- -mean(PCAWG_target[loc + length(list_targets_chrom)])
    # }
    # #   Add fitted parameters to Davoli-Charm table
    # CHARM_table$Selection_rate <- 0
    # for (row in 1:nrow(CHARM_table)) {
    #     arm <- CHARM_table$Arm[row]
    #     loc <- which(list_parameters_arm == arm)
    #     CHARM_table$Selection_rate[row] <- parameters_best[loc]
    # }
    # #   Add fitted parameters to Davoli-Chrom table
    # CHROM_table$Selection_rate <- 0
    # for (row in 1:nrow(CHROM_table)) {
    #     chrom <- CHROM_table$Chromosome[row]
    #     loc <- which(list_parameters_chrom == chrom)
    #     CHROM_table$Selection_rate[row] <- mean(parameters_best[loc])
    # }
    # #   Plot comparison between fitted arm-level selection rates and Davoli Charm scores
    # filename <- paste(model_name, "_against_Davoli.jpeg", sep = "")
    # plot_Charm(CHARM_table, CHROM_table, filename, model_name)
}

# plot_Charm <- function(CHARM_table, CHROM_table, filename, model_name) {
#     print(CHARM_table)
#     print(CHROM_table)
#     #-----------------------------------------Configuration for subplots
#     layout <- matrix(NA, nrow = 3, ncol = 2)
#     gs <- list()
#     plots <- list()
#     id <- 0
#     #------------Plot selection rates vs cancer-specific amp frequencies
#     id <- id + 1
#     layout[1, 1] <- id
#     tmp <- cor.test(CHARM_table$Selection_rate, CHARM_table$Amp_freq_spec)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.1, y = 0.9, just = "left",
#             gp = gpar(col = "blue", fontsize = 30, fontface = "bold"),
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Selection_rate, y = Amp_freq_spec)) +
#         geom_point(size = 25, color = "blue") +
#         geom_smooth(method = lm, color = "blue") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Selection rates") +
#         ylab(paste("Amplification frequencies (", model_name, ")", sep = "")) +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #------------Plot selection rates vs cancer-specific del frequencies
#     id <- id + 1
#     layout[1, 2] <- id
#     tmp <- cor.test(CHARM_table$Selection_rate, CHARM_table$Del_freq_spec)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.9, y = 0.9, just = "right",
#             gp = gpar(col = "blue", fontsize = 30, fontface = "bold")
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Selection_rate, y = Del_freq_spec)) +
#         geom_point(size = 25, color = "blue") +
#         geom_smooth(method = lm, color = "blue") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Selection rates") +
#         ylab(paste("Deletion frequencies (", model_name, ")", sep = "")) +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #-------------------------Plot selection rates vs Charm.TSG.OG score
#     id <- id + 1
#     layout[2, 1] <- id
#     tmp <- cor.test(CHARM_table$Selection_rate, CHARM_table$Charm.TSG.OG.score)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.1, y = 0.9, just = "left",
#             gp = gpar(col = "coral", fontsize = 30, fontface = "bold")
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Selection_rate, y = Charm.TSG.OG.score)) +
#         geom_point(size = 25, color = "coral") +
#         geom_smooth(method = lm, color = "coral") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Selection rates") +
#         ylab("Charm(TSG,OG) score") +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #---------------------Plot selection rates vs Charm.TSG.OG.Ess score
#     id <- id + 1
#     layout[2, 2] <- id
#     tmp <- cor.test(CHARM_table$Selection_rate, CHARM_table$Charm.TSG.OG.Ess.score)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.9, y = 0.1, just = "right",
#             gp = gpar(col = "coral", fontsize = 30, fontface = "bold")
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Selection_rate, y = Charm.TSG.OG.Ess.score)) +
#         geom_point(size = 25, color = "coral") +
#         geom_smooth(method = lm, color = "coral") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Selection rates") +
#         ylab("Charm(TSG,OG,Ess) score") +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #-------------Plot cancer-specific amp freqs vs pan-cancer amp freqs
#     id <- id + 1
#     layout[3, 1] <- id
#     tmp <- cor.test(CHARM_table$Amp_freq_all, CHARM_table$Amp_freq_spec)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.9, y = 0.1, just = "right",
#             gp = gpar(col = "darkgreen", fontsize = 30, fontface = "bold")
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Amp_freq_all, y = Amp_freq_spec)) +
#         geom_point(size = 25, color = "darkgreen") +
#         geom_smooth(method = lm, color = "darkgreen") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Amplification frequencies (pan-cancer)") +
#         ylab(paste("Amplification frequencies (", model_name, ")", sep = "")) +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #-------------Plot cancer-specific del freqs vs pan-cancer del freqs
#     id <- id + 1
#     layout[3, 2] <- id
#     tmp <- cor.test(CHARM_table$Del_freq_all, CHARM_table$Del_freq_spec)
#     Pearson_r <- tmp$estimate
#     Pearson_p_val <- tmp$p.value
#     grob1 <- grobTree(
#         textGrob(
#             paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
#             x = 0.9, y = 0.1, just = "right",
#             gp = gpar(col = "darkgreen", fontsize = 30, fontface = "bold")
#         )
#     )
#     plots[[id]] <- ggplot(CHARM_table, aes(x = Del_freq_all, y = Del_freq_spec)) +
#         geom_point(size = 25, color = "darkgreen") +
#         geom_smooth(method = lm, color = "darkgreen") +
#         annotation_custom(grob1) +
#         geom_text(aes(label = Arm), size = 10, color = "white") +
#         xlab("Deletion frequencies (pan-cancer)") +
#         ylab(paste("Deletion frequencies (", model_name, ")", sep = "")) +
#         # xlim(c(0.5, 1.5)) +
#         # ylim(c(0, 1)) +
#         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#         theme(text = element_text(size = 30))
#     gs[[id]] <- plots[[id]]
#     #----------------------------------------------------Print all plots
#     jpeg(filename, width = 3000, height = 2000)
#     p <- grid.arrange(grobs = gs, layout_matrix = layout)
#     print(p)
#     dev.off()
# }
