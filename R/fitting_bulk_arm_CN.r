#----------------------Function to assign parameters to proper positions
assign_paras <- function(model_variables, parameter_IDs, parameters) {
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
#--------------------------Function to extract arm-level gain/loss delta
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
#-------------------------------------Objective function for ABC fitting
func_ABC <- function(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets) {
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile"),
        R_libPaths = R_libPaths
    )
    #   Statistics = arm-level gain/loss delta
    SIMS_delta_genome_arms <- gainloss_SIMS(SIMS_chromosome, ploidy_normalization = FALSE)
    SIMS_arm_gainloss <- get_arm_gainloss(SIMS_delta_genome_arms, copynumber_coordinates, list_targets)
    stat <- c(SIMS_arm_gainloss$delta_gain, SIMS_arm_gainloss$delta_loss)
}
#-----------------Function to choose one best parameter from a posterior
get_best_para <- function(data_rf, model_rf, obs_rf, post_rf) {
    df_dist <- densityPlot_df(model_rf, obs_rf, data_rf)
    best_para <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
}

#' @export
library_bulk_arm_CN <- function(library_name,
                                model_variables,
                                list_parameters,
                                list_targets,
                                ABC_simcount = 10000,
                                n_cores = NULL,
                                n_samples = 100,
                                R_libPaths = NULL) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    #---------------------------------List of parameter IDs to be fitted
    parameter_IDs <- list_parameters$Variable
    #---------------------Find arm-level gain/loss map for entire genome
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
    # =============================================CREATE REFERENCE DATA
    #---------------------------------------Simulate table of parameters
    sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    for (col in 1:ncol(sim_param)) {
        sim_param[, col] <- runif(ABC_simcount, min = as.numeric(list_parameters$Lower_bound[col]), max = as.numeric(list_parameters$Upper_bound[col]))
    }
    #-----------------------------------------------Make reference table
    start_time <- Sys.time()
    #   Configure parallel pool
    cl <- makePSOCKcluster(n_cores)
    cat("Creating reference table for ABC...\n")
    n_samples <<- n_samples
    list_targets <<- list_targets
    sim_param <<- sim_param
    parameter_IDs <<- parameter_IDs
    model_variables <<- model_variables
    gainloss_SIMS <<- gainloss_SIMS
    func_ABC <<- func_ABC
    assign_paras <<- assign_paras
    get_arm_gainloss <<- get_arm_gainloss
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c(
        "n_samples", "list_targets", "sim_param", "parameter_IDs", "model_variables", "gainloss_SIMS", "func_ABC", "assign_paras", "get_arm_gainloss", "R_libPaths",
        "BUILD_driver_library", "simulator_full_program", "one_simulation",
        "SIMULATOR_VARIABLES_for_simulation",
        "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
        "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
        "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
        "SIMULATOR_FULL_PHASE_2_main",
        "SIMULATOR_FULL_PHASE_3_main",
        "get_cn_profile", "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model"
    ))
    e <- new.env()
    e$libs <- .libPaths()
    clusterExport(cl, "libs", envir = e)
    clusterEvalQ(cl, .libPaths(libs))
    #   Create simulated statistics in parallel
    pbo <- pboptions(type = "txt")
    sim_results_list <- pblapply(cl = cl, X = 1:ABC_simcount, FUN = function(iteration) {
        parameters <- sim_param[iteration, ]
        stat <- func_ABC(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets)
        return(stat)
    })
    stopCluster(cl)
    #   Group simulated statistics into one table
    sim_stat <- matrix(0, nrow = ABC_simcount, ncol = 2 * length(list_targets))
    for (row in 1:ABC_simcount) {
        stat <- sim_results_list[[row]]
        sim_stat[row, ] <- stat
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    #---------------------------Save the parameters and their statistics
    ABC_input <- list()
    ABC_input$model_variables <- model_variables
    ABC_input$n_samples <- n_samples
    ABC_input$parameter_IDs <- parameter_IDs
    ABC_input$sim_param <- sim_param
    ABC_input$sim_stat <- sim_stat
    filename <- paste0(library_name, "_ABC_input.rda")
    save(ABC_input, file = filename)
}

#' @export
fitting_bulk_arm_CN <- function(library_name,
                                model_name,
                                copynumber_DATA,
                                type_sample_DATA = "individual",
                                type_cn_DATA = "bin",
                                list_parameters,
                                list_targets,
                                n_cores = NULL,
                                R_libPaths = NULL,
                                folder_workplace = NULL) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    folder_workplace_tmp <- folder_workplace
    if (is.null(folder_workplace_tmp)) {
        folder_workplace_tmp <- ""
    } else {
        dir.create(folder_workplace_tmp)
        folder_workplace_tmp <- paste(folder_workplace_tmp, "/", sep = "")
    }
    #-----------------------------------------Input simulated CN library
    filename <- paste0(library_name, "_ABC_input.rda")
    load(filename)
    model_variables <- ABC_input$model_variables
    n_samples <- ABC_input$n_samples
    parameter_IDs <- ABC_input$parameter_IDs
    sim_param <- ABC_input$sim_param
    sim_stat <- ABC_input$sim_stat
    #---------------------Find arm-level gain/loss map for entire genome
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
    #   Get the target statistics (= gain/loss on each chromosome arm)
    if (type_sample_DATA == "individual") {
        #   Get gain/loss map for this chromosome at arm level
        DATA_delta_genome_arms <- gainloss_DATA(
            copynumber_DATA, copynumber_coordinates,
            ploidy_normalization = TRUE,
            use_rbindlist = TRUE,
            arm_level = TRUE, pos_centromeres = cn_info
        )
        #   Get the target statistics (= gain/loss on each chromosome arm)
        DATA_arm_gainloss <- get_arm_gainloss(DATA_delta_genome_arms, copynumber_coordinates, list_targets)
        DATA_target <- c(DATA_arm_gainloss$delta_gain, DATA_arm_gainloss$delta_loss)
    } else if (type_sample_DATA == "average") {
        DATA_target <- rep(0, 2 * length(list_targets))
        for (i in 1:length(list_targets)) {
            arm <- list_targets[i]
            loc <- which(copynumber_DATA$Arm == arm)
            if (length(loc) > 0) {
                DATA_target[i] <- copynumber_DATA$Amp_freq_all[loc]
                DATA_target[i + length(list_targets)] <- -copynumber_DATA$Del_freq_all[loc]
            }
        }
    }
    # # ====================================FITTING WITH ABC RANDOM FOREST
    # #--------------------------------------------Fit parameters with ABC
    # #---Dataframe for data observation
    # obs_rf <- data.frame(matrix(DATA_target, nrow = 1))
    # colnames(obs_rf) <- paste0("gainloss_", 1:ncol(obs_rf))
    # #---Dataframe for parameters for reference
    # all_paras <- data.frame(sim_param[, which(parameter_IDs %in% list_parameters$Variable)])
    # colnames(all_paras) <- parameter_IDs[which(parameter_IDs %in% list_parameters$Variable)]
    # #---Dataframe for corresponding statistics for reference
    # all_data <- data.frame(sim_stat[
    #     , c(
    #         which(model_variables$chromosome_arm_library$Arm_ID %in% list_targets),
    #         which(model_variables$chromosome_arm_library$Arm_ID %in% list_targets) + length(model_variables$chromosome_arm_library$Arm_ID)
    #     )
    # ])
    # colnames(all_data) <- paste0("gainloss_", 1:ncol(all_data))
    # #---Fit each parameter with ABC-rf
    # layout <- matrix(NA, nrow = 7, ncol = ceiling(length(parameter_IDs) / 7))
    # gs <- list()
    # id <- 0
    # # for (para in 1:length(parameter_IDs)) {
    # for (para in 1:nrow(list_parameters)) {
    #     para_ID <- list_parameters$Variable[para]
    #     # para_ID <- parameter_IDs[para]
    #     cat(paste("ABC for parameter ", para_ID, " [", para, "/", nrow(list_parameters), "]", "\n", sep = ""))
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
    #     filename <- paste0(folder_workplace_tmp, model_name, "_ABC_output_", para_ID, ".rda")
    #     save(ABC_output, file = filename)
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
    list_parameters$Best_value <- 0
    for (para in 1:nrow(list_parameters)) {
        para_ID <- list_parameters$Variable[para]
        filename <- paste0(folder_workplace_tmp, model_name, "_ABC_output_", para_ID, ".rda")
        load(filename)
        best_rf <- ABC_output$best_rf
        parameters_best[para] <- best_rf
        list_parameters$Best_value[para] <- best_rf
        cat(paste0(para_ID, "===", best_rf, "\n"))
    }
    filename <- paste0(model_name, "_fitted_parameters.csv")
    write.csv(list_parameters, filename)



    #------------------------Analysis of fitted CN profiles against data
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, list_parameters$Variable, parameters_best)
    print(model_variables$chromosome_arm_library)
    # model_variables <- assign_paras(model_variables, parameter_IDs, parameters_best)
    #   Make simulations using best parameters
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE, compute_parallel = TRUE, n_cores = n_cores,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile"),
        R_libPaths = R_libPaths
    )
    #   Plot comparison vs bin-level data
    if (type_cn_DATA == "bin") {
        filename <- paste(model_name, "_against_bin_data.jpeg", sep = "")
        plot_gainloss(
            copynumber_sims = SIMS_chromosome, copynumber_DATA = copynumber_DATA,
            title = model_name, filename = filename,
            type_sample_DATA = type_sample_DATA,
            arm_level = FALSE, pos_centromeres = cn_info
        )
    }
    #   Plot comparison vs arm-level data
    filename <- paste(model_name, "_against_arm_data.jpeg", sep = "")
    plot_gainloss(
        copynumber_sims = SIMS_chromosome, copynumber_DATA = copynumber_DATA,
        title = model_name, filename = filename,
        type_sample_DATA = type_sample_DATA,
        arm_level = TRUE, pos_centromeres = cn_info
    )
    #----------------Analysis of fitted parameters against amp/del freqs
    #   Prepare plot dataframe with fitted selection rates and amp/del freqs
    plot_table <- data.frame(Arm = list_targets)
    plot_table$Amp_freq_spec <- 0
    plot_table$Del_freq_spec <- 0
    plot_table$Selection_rate <- 0
    for (row in 1:nrow(plot_table)) {
        arm <- plot_table$Arm[row]
        plot_table$Amp_freq_spec[row] <- DATA_target[which(list_targets == arm)]
        plot_table$Del_freq_spec[row] <- -DATA_target[which(list_targets == arm) + length(list_targets)]
        plot_table$Selection_rate[row] <- parameters_best[which(list_parameters$Variable == arm)]
    }
    #   Configuration for subplots
    layout <- matrix(NA, nrow = 1, ncol = 2)
    gs <- list()
    plots <- list()
    id <- 0
    #   Plot selection rates vs cancer-specific amp frequencies
    id <- id + 1
    layout[1, 1] <- id
    tmp <- cor.test(plot_table$Selection_rate, plot_table$Amp_freq_spec)
    Pearson_r <- tmp$estimate
    Pearson_p_val <- tmp$p.value
    grob1 <- grobTree(
        textGrob(
            paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
            x = 0.1, y = 0.9, just = "left",
            gp = gpar(col = "blue", fontsize = 50, fontface = "bold"),
        )
    )
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Amp_freq_spec)) +
        geom_point(size = 30, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        annotation_custom(grob1) +
        geom_text(aes(label = Arm), size = 12, color = "white") +
        xlab("Selection rates") +
        ylab("Amplification frequencies") +
        # ylab(paste("Amplification frequencies (", model_name, ")", sep = "")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 70)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    gs[[id]] <- plots[[id]]
    #   Plot selection rates vs cancer-specific del frequencies
    id <- id + 1
    layout[1, 2] <- id
    tmp <- cor.test(plot_table$Selection_rate, plot_table$Del_freq_spec)
    Pearson_r <- tmp$estimate
    Pearson_p_val <- tmp$p.value
    grob1 <- grobTree(
        textGrob(
            paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
            x = 0.9, y = 0.9, just = "right",
            gp = gpar(col = "blue", fontsize = 50, fontface = "bold")
        )
    )
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Del_freq_spec)) +
        geom_point(size = 30, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        annotation_custom(grob1) +
        geom_text(aes(label = Arm), size = 12, color = "white") +
        xlab("Selection rates") +
        ylab("Deletion frequencies") +
        # ylab(paste("Deletion frequencies (", model_name, ")", sep = "")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 70)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    gs[[id]] <- plots[[id]]
    #   Print all plots
    filename <- paste(model_name, "_selection_rate_analysis.jpeg", sep = "")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
    #----------------Analysis of fitted parameters against Davoli scores
    if ("Charm.TSG.OG.Ess.score" %in% colnames(copynumber_DATA)) {
        #   Supplement plot dataframe with Charm scores
        plot_table$Charm.TSG.OG.score <- NA
        plot_table$Charm.TSG.OG.Ess.score <- NA
        for (i in 1:nrow(copynumber_DATA)) {
            arm <- copynumber_DATA$Arm[i]
            loc <- which(plot_table$Arm == arm)
            plot_table$Charm.TSG.OG.score[loc] <- copynumber_DATA$Charm.TSG.OG.score[i]
            plot_table$Charm.TSG.OG.Ess.score[loc] <- copynumber_DATA$Charm.TSG.OG.Ess.score[i]
        }
        if (any(is.na(plot_table$Charm.TSG.OG.Ess.score))) plot_table <- plot_table[-which(is.na(plot_table$Charm.TSG.OG.Ess.score)), ]
        #   Configuration for subplots
        layout <- matrix(NA, nrow = 1, ncol = 2)
        gs <- list()
        plots <- list()
        id <- 0
        #   Plot selection rates vs Charm.TSG.OG score
        id <- id + 1
        layout[1, 1] <- id
        tmp <- cor.test(plot_table$Selection_rate, plot_table$Charm.TSG.OG.score)
        Pearson_r <- tmp$estimate
        Pearson_p_val <- tmp$p.value
        grob1 <- grobTree(
            textGrob(
                paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
                x = 0.9, y = 0.9, just = "right",
                gp = gpar(col = "coral", fontsize = 50, fontface = "bold")
            )
        )
        plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Charm.TSG.OG.score)) +
            geom_point(size = 30, color = "coral") +
            geom_smooth(method = lm, color = "coral") +
            annotation_custom(grob1) +
            geom_text(aes(label = Arm), size = 12, color = "white") +
            xlab("Selection rates") +
            ylab("Charm(TSG,OG) score") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 70)) +
            theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
        gs[[id]] <- plots[[id]]
        #   Plot selection rates vs Charm.TSG.OG.Ess score
        id <- id + 1
        layout[1, 2] <- id
        tmp <- cor.test(plot_table$Selection_rate, plot_table$Charm.TSG.OG.Ess.score)
        Pearson_r <- tmp$estimate
        Pearson_p_val <- tmp$p.value
        grob1 <- grobTree(
            textGrob(
                paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
                x = 0.9, y = 0.9, just = "right",
                gp = gpar(col = "coral", fontsize = 50, fontface = "bold")
            )
        )
        plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Charm.TSG.OG.Ess.score)) +
            geom_point(size = 30, color = "coral") +
            geom_smooth(method = lm, color = "coral") +
            annotation_custom(grob1) +
            geom_text(aes(label = Arm), size = 12, color = "white") +
            xlab("Selection rates") +
            ylab("Charm(TSG,OG,Ess) score") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 70)) +
            theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
        gs[[id]] <- plots[[id]]
        #   Print all plots
        filename <- paste(model_name, "_selection_rate_vs_Charm_scores.jpeg", sep = "")
        jpeg(filename, width = 3000, height = 1500)
        p <- grid.arrange(grobs = gs, layout_matrix = layout)
        print(p)
        dev.off()
    }
}

#' @export
statistics_bulk_arm_WGD_against_losses <- function(plotname,
                                                   DATA_cancer_types,
                                                   DATA_cancer_type_sample_ids,
                                                   DATA_wgd) {
    library(ggplot2)
    library(ggrepel)

    #---------------------------Find WGD proportions in each cancer type
    DATA_wgd_proportion <- rep(0, length(DATA_cancer_types))
    for (i in 1:length(DATA_cancer_types)) {
        mini_DATA_wgd <- DATA_wgd[which(DATA_wgd$samplename %in% DATA_cancer_type_sample_ids[[i]] & DATA_wgd$wgd_uncertain == FALSE), ]
        DATA_wgd_proportion[i] <- length(which(mini_DATA_wgd$wgd_status == "wgd")) / length(mini_DATA_wgd$wgd_status)
    }
    #------------Find count and strength of TSG arms in each cancer type
    FIT_tsg_count <- rep(0, length(DATA_cancer_types))
    FIT_tsg_strength <- rep(0, length(DATA_cancer_types))
    for (i in 1:length(DATA_cancer_types)) {
        cancer_types_fit <- read.csv(paste0(DATA_cancer_types[i], "_fitted_parameters.csv"), header = TRUE)
        cancer_types_fit_selection_rates <- cancer_types_fit$Best_value[which(cancer_types_fit$Type == "Arm_selection_rate")]
        FIT_tsg_count[i] <- length(which(cancer_types_fit_selection_rates < 1))
        if (FIT_tsg_count[i] == 0) {
            FIT_tsg_strength[i] <- 1
        } else {
            FIT_tsg_strength[i] <- 1 / mean(cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates < 1)])
        }
    }
    #-------------------Plot relationship between fitted selection rates
    #---------------------------------------and observed WGD proportions
    filename <- paste0(plotname, "_WGD_vs_fitted_TSG.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    df_plot <- data.frame(
        cancer_types = DATA_cancer_types,
        WGD = DATA_wgd_proportion,
        TSG_count = FIT_tsg_count,
        TSG_strength = FIT_tsg_strength
    )
    p <- ggplot(df_plot, aes(x = TSG_count, y = TSG_strength, color = WGD)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_types), size = 10, point.padding = 0.4) +
        xlab("Count of TSG arms") +
        ylab("Mean selection rate of TSG arms") +
        labs(fill = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()


    tmp <- cor.test(FIT_tsg_count, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    print("COUNT:")
    print(tmp$estimate)
    print(tmp$p.value)

    tmp <- cor.test(FIT_tsg_strength, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    print("STRENGTH:")
    print(tmp$estimate)
    print(tmp$p.value)


    # print(DATA_cancer_type_sample_ids)
    # print(DATA_wgd)
    print("-----------------------------------------------------------")
    print(DATA_cancer_types)
    print("-----------------------------------------------------------")
    print(DATA_wgd_proportion)
    print("-----------------------------------------------------------")
    print(FIT_tsg_count)
    print("-----------------------------------------------------------")
    print(FIT_tsg_strength)
}
