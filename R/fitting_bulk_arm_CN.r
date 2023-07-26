#----------------------Function to assign parameters to proper positions
#' @export
bulk_arm_CN_assign_paras <- function(model_variables, parameter_IDs, parameters) {
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
bulk_arm_CN_get_arm_gainloss <- function(bin_gainloss, copynumber_coordinates, list_targets) {
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
#-----Objective function for ABC fitting for CNA rates & selection rates
bulk_arm_CN_func_ABC <- function(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets) {
    #   Assign parameters in model variables
    model_variables <- bulk_arm_CN_assign_paras(model_variables, parameter_IDs, parameters)
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile"),
        R_libPaths = R_libPaths
    )
    #   Statistics = arm-level gain/loss delta
    SIMS_delta_genome_arms <- gainloss_SIMS(SIMS_chromosome, ploidy_normalization = FALSE)
    SIMS_arm_gainloss <- bulk_arm_CN_get_arm_gainloss(SIMS_delta_genome_arms, copynumber_coordinates, list_targets)
    stat <- c(SIMS_arm_gainloss$delta_gain, SIMS_arm_gainloss$delta_loss)
}
#-----------------Function to choose one best parameter from a posterior
bulk_arm_CN_get_best_para <- function(data_rf, model_rf, obs_rf, post_rf) {
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
    # gainloss_SIMS <<- gainloss_SIMS
    bulk_arm_CN_func_ABC <<- bulk_arm_CN_func_ABC
    bulk_arm_CN_assign_paras <<- bulk_arm_CN_assign_paras
    bulk_arm_CN_get_arm_gainloss <<- bulk_arm_CN_get_arm_gainloss
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c(
        "n_samples", "list_targets", "sim_param", "parameter_IDs", "model_variables", "gainloss_SIMS", "bulk_arm_CN_func_ABC", "bulk_arm_CN_assign_paras", "bulk_arm_CN_get_arm_gainloss", "R_libPaths",
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
        stat <- bulk_arm_CN_func_ABC(parameters, parameter_IDs, model_variables, copynumber_coordinates, list_targets)
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
                                model_variables,
                                copynumber_DATA,
                                type_sample_DATA = "individual",
                                type_cn_DATA = "bin",
                                list_parameters,
                                list_parameters_library,
                                list_targets,
                                list_targets_library,
                                bound_freq = 0.1,
                                ntree = 200,
                                library_shuffle = FALSE,
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
    library(data.table)
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
    # model_variables <- ABC_input$model_variables
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
        DATA_arm_gainloss <- bulk_arm_CN_get_arm_gainloss(DATA_delta_genome_arms, copynumber_coordinates, list_targets)
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
    # ==========================INCREASE SIMULATED LIBRARY VIA BOOTSTRAP
    #   Find ID for each parameter in the prepared library
    sim_param_ID <- list_parameters_library$Variable
    sim_param_type <- list_parameters_library$Type
    sim_stat_ID <- c(paste0("gain_", list_targets_library), paste0("loss_", list_targets_library))
    #   Find unique chromosomes in the library
    sim_stat_chroms <- sim_stat_ID
    for (i in 1:length(sim_stat_chroms)) {
        sim_stat_chroms[i] <- sub(".*_", "", substr(sim_stat_chroms[i], 1, nchar(sim_stat_chroms[i]) - 1))
    }
    unique_chroms <- unique(sim_stat_chroms)
    #   Prepare the original prepared library
    df_sim_param <- data.frame(sim_param)
    colnames(df_sim_param) <- sim_param_ID
    df_sim_stat <- data.frame(sim_stat)
    colnames(df_sim_stat) <- sim_stat_ID
    #   Bootstrap each simulation in the prepared library
    ls_sim_param_bootstrap <- vector("list", 1)
    ls_sim_param_bootstrap[[1]] <- df_sim_param
    ls_sim_stat_bootstrap <- vector("list", 1)
    ls_sim_stat_bootstrap[[1]] <- df_sim_stat
    tmp <- 0
    if (library_shuffle) {
        N_shuffle <- length(unique_chroms) - 1
        for (i in 1:N_shuffle) {
            df_sim_param_next <- df_sim_param
            df_sim_stat_next <- df_sim_stat
            #   Shuffle the chromosome indices
            shuffle_chroms <- c(unique_chroms[(i + 1):length(unique_chroms)], unique_chroms[1:i])
            #   Shuffle the library of parameters accordingly
            for (j in 1:ncol(df_sim_param)) {
                if (sim_param_type[j] != "Arm_selection_rate") next
                param_ID_old <- sim_param_ID[j]
                param_chrom_old <- substr(param_ID_old, 1, nchar(param_ID_old) - 1)
                param_chrom_new <- shuffle_chroms[which(unique_chroms == param_chrom_old)]
                param_ID_new <- paste0(param_chrom_new, substr(param_ID_old, nchar(param_ID_old), nchar(param_ID_old)))
                df_sim_param_next[param_ID_new] <- df_sim_param[param_ID_old]
            }
            #   Shuffle the library of statistics accordingly
            for (j in 1:ncol(df_sim_stat)) {
                param_ID_old <- sim_stat_ID[j]
                param_chrom_old <- sub(".*_", "", substr(param_ID_old, 1, nchar(param_ID_old) - 1))
                param_chrom_new <- shuffle_chroms[which(unique_chroms == param_chrom_old)]
                param_ID_new <- paste0(sub("\\_.*", "_", param_ID_old), param_chrom_new, substr(param_ID_old, nchar(param_ID_old), nchar(param_ID_old)))
                df_sim_stat_next[param_ID_new] <- df_sim_stat[param_ID_old]
            }
            #   Record the shuffled library
            ls_sim_param_bootstrap[[i + 1]] <- df_sim_param_next
            ls_sim_stat_bootstrap[[i + 1]] <- df_sim_stat_next
        }
    }
    #   Combine all shuffled libraries into one dataframe
    sim_param_bootstrap <- as.data.frame(rbindlist(ls_sim_param_bootstrap))
    sim_stat_bootstrap <- as.data.frame(rbindlist(ls_sim_stat_bootstrap))
    # ====================================FITTING WITH ABC RANDOM FOREST
    #---Dataframe for prepared library of parameters
    all_paras_original <- df_sim_param
    all_paras_bootstrap <- sim_param_bootstrap
    #---Dataframe for prepared library of statistics
    all_data_original <- df_sim_stat
    all_data_bootstrap <- sim_stat_bootstrap
    #---Dataframe for data observation
    all_obs <- data.frame(matrix(DATA_target, nrow = 1))
    colnames(all_obs) <- c(paste0("gain_", list_targets), paste0("loss_", list_targets))
    #---Fit each parameter with ABC-rf
    layout <- matrix(NA, nrow = 7, ncol = ceiling(length(parameter_IDs) / 7))
    gs <- list()
    id <- 0
    for (para in 1:nrow(list_parameters)) {
        start_time <- Sys.time()
        para_ID <- list_parameters$Variable[para]
        para_type <- list_parameters$Type[para]
        cat(paste("\nABC for parameter ", para_ID, " [", para, "/", nrow(list_parameters), "]", "\n", sep = ""))
        #   Prepare observations for this parameter
        if (para_type == "CNA_probability") {
            mini_obs <- all_obs
        } else if (para_type == "Arm_selection_rate") {
            para_arm <- para_ID
            para_chrom <- substr(para_ID, 1, nchar(para_ID) - 1)
            mini_obs <- NULL
            for (stat in colnames(all_obs)) {
                stat_arm <- sub(".*_", "", stat)
                stat_chrom <- substr(stat_arm, 1, nchar(stat_arm) - 1)
                if (stat_chrom != para_chrom) next
                if (is.null(mini_obs)) {
                    mini_obs <- all_obs[stat]
                } else {
                    mini_obs <- cbind(mini_obs, all_obs[stat])
                }
            }
        }
        #   Prepare library of statistics for this parameter
        mini_data <- NULL
        for (stat in colnames(mini_obs)) {
            if (para_type == "CNA_probability") {
                next_data <- all_data_original[stat]
            } else if (para_type == "Arm_selection_rate") {
                next_data <- all_data_bootstrap[stat]
            }
            if (is.null(mini_data)) {
                mini_data <- next_data
            } else {
                mini_data <- cbind(mini_data, next_data)
            }
        }
        #   Prepare library of parameters for this parameter
        if (para_type == "CNA_probability") {
            data_rf <- cbind(all_paras_original[para_ID], mini_data)
        } else if (para_type == "Arm_selection_rate") {
            data_rf <- cbind(all_paras_bootstrap[para_ID], mini_data)
        }
        #   No fitting for chromosomes without frequent gains/losses
        if ((bound_freq > 0) & (para_type == "Arm_selection_rate")) {
            #   Find condition of data
            para_all_arms <- unique(sub(".*_", "", colnames(mini_obs)))
            para_all_arms_classification <- rep(0, length(para_all_arms))
            for (i in 1:length(para_all_arms)) {
                if (mini_obs[[paste0("gain_", para_all_arms[i])]] >= -mini_obs[[paste0("loss_", para_all_arms[i])]]) {
                    para_all_arms_classification[i] <- mini_obs[[paste0("gain_", para_all_arms[i])]]
                } else {
                    para_all_arms_classification[i] <- mini_obs[[paste0("loss_", para_all_arms[i])]]
                }
            }
            #   If both arms are not significantly gains/lost, then no fitting
            if (max(abs(para_all_arms_classification) < bound_freq)) {
                best_rf <- 1
                #   Save results for fitting this parameter
                ABC_output <- list()
                ABC_output$para_ID <- para_ID
                ABC_output$best_rf <- best_rf
                filename <- paste0(folder_workplace_tmp, model_name, "_ABC_output_", para_ID, ".rda")
                save(ABC_output, file = filename)
                cat(paste0("Best parameter: ", best_rf, "\n"))
                next
            }
        }
        #   Train the random forest
        colnames(data_rf)[1] <- "para"
        f <- as.formula("para ~.")
        model_rf <- regAbcrf(
            formula = f, data_rf,
            paral = TRUE, ncores = n_cores,
            ntree = ntree,
            sampsize = nrow(data_rf),
            save.memory = TRUE
        )
        #   Predict posterior distribution based on found random forest
        post_rf <- predict(model_rf, mini_obs, data_rf, paral = TRUE, ncores = n_cores)
        #   Choose best value from posterior distribution
        best_rf <- bulk_arm_CN_get_best_para(data_rf, model_rf, mini_obs, post_rf)
        #   Save results for fitting this parameter
        ABC_output <- list()
        ABC_output$para_ID <- para_ID
        ABC_output$best_rf <- best_rf
        filename <- paste0(folder_workplace_tmp, model_name, "_ABC_output_", para_ID, ".rda")
        save(ABC_output, file = filename)
        #   Plot the prior, posterior and chosen best parameter for all variables
        id <- id + 1
        row <- id %% 7
        if (row == 0) row <- 7
        col <- ceiling(id / 7)
        layout[row, col] <- id
        gs[[id]] <- densityPlot_MODIFIED(
            model_rf, all_obs, data_rf,
            protocol = "arm",
            fontsize = 20,
            chosen_para = best_rf,
            color_prior = "lightblue", color_posterior = "darkblue", color_vline = "blue",
            main = para_ID
        )
        end_time <- Sys.time()
        cat(paste0("Best parameter: ", best_rf, "\n"))
        print(end_time - start_time)
        #   Clear memory
        model_rf <- c()
        data_rf <- c()
        post_rf <- c()
        best_rf <- c()
    }
    #   Plot the prior, posterior and chosen best parameter for all variables
    filename <- paste0(model_name, "_ABC_all.jpeg")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
    # =======================================ANALYSIS OF FITTING RESULTS
    #------------------Choose the best parameter set from all posteriors
    parameter_IDs_best <- list_parameters$Variable
    parameters_best <- rep(0, nrow(list_parameters))
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
    model_variables <- bulk_arm_CN_assign_paras(model_variables, parameter_IDs_best, parameters_best)
    print(model_variables$chromosome_arm_library)
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
        plot_table$Selection_rate[row] <- parameters_best[which(parameter_IDs_best == arm)]
    }
    #   Remove chromosome arms that were not fitted
    if (length(which(plot_table$Selection_rate == 1)) > 0) {
        plot_table <- plot_table[-which(plot_table$Selection_rate == 1), ]
    }
    #   Configuration for subplots
    layout <- matrix(NA, nrow = 1, ncol = 2)
    gs <- list()
    plots <- list()
    id <- 0
    #   Plot selection rates vs cancer-specific amp frequencies
    id <- id + 1
    layout[1, 1] <- id
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Amp_freq_spec)) +
        geom_point(size = 30, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        geom_text(aes(label = Arm), size = 12, color = "white") +
        xlab("Selection rates") +
        ylab("Amplification frequencies") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 70)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    if (nrow(plot_table) >= 5) {
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
        plots[[id]] <- plots[[id]] + annotation_custom(grob1)
    }
    gs[[id]] <- plots[[id]]
    #   Plot selection rates vs cancer-specific del frequencies
    id <- id + 1
    layout[1, 2] <- id
    plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Del_freq_spec)) +
        geom_point(size = 30, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        geom_text(aes(label = Arm), size = 12, color = "white") +
        xlab("Selection rates") +
        ylab("Deletion frequencies") +
        # ylab(paste("Deletion frequencies (", model_name, ")", sep = "")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 70)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    if (nrow(plot_table) >= 5) {
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
        plots[[id]] <- plots[[id]] + annotation_custom(grob1)
    }
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
        plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Charm.TSG.OG.score)) +
            geom_point(size = 30, color = "coral") +
            geom_smooth(method = lm, color = "coral") +
            geom_text(aes(label = Arm), size = 12, color = "white") +
            xlab("Selection rates") +
            ylab("Charm(TSG,OG) score") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 70)) +
            theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
        if (nrow(plot_table) >= 5) {
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
            plots[[id]] <- plots[[id]] + annotation_custom(grob1)
        }
        gs[[id]] <- plots[[id]]
        #   Plot selection rates vs Charm.TSG.OG.Ess score
        id <- id + 1
        layout[1, 2] <- id
        plots[[id]] <- ggplot(plot_table, aes(x = Selection_rate, y = Charm.TSG.OG.Ess.score)) +
            geom_point(size = 30, color = "coral") +
            geom_smooth(method = lm, color = "coral") +
            geom_text(aes(label = Arm), size = 12, color = "white") +
            xlab("Selection rates") +
            ylab("Charm(TSG,OG,Ess) score") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 70)) +
            theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
        if (nrow(plot_table) >= 5) {
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
            plots[[id]] <- plots[[id]] + annotation_custom(grob1)
        }
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
statistics_bulk_arm_WGD_status <- function(plotname,
                                           DATA_cancer_types,
                                           DATA_cancer_type_sample_ids,
                                           DATA_cancer_type_cn,
                                           DATA_wgd,
                                           model_variables) {
    library(ggplot2)
    library(ggrepel)
    library(scales)
    #-------------------------Find genome coordinate from one simulation
    cn_info <- model_variables$cn_info
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
    #----------------Find WGD rate & FGA difference for each cancer type
    df_WGD_FGA <- data.frame(matrix(0, nrow = length(DATA_cancer_types), ncol = 3))
    colnames(df_WGD_FGA) <- c("cancer_type", "WGD", "WGD_increased_FGA")
    pb <- txtProgressBar(min = 0, max = length(DATA_cancer_types), style = 3, width = 50, char = "=")
    # for (i in 1:length(DATA_cancer_types)) {
    for (i in 1:2) {
        setTxtProgressBar(pb, i)
        cancer_type <- DATA_cancer_types[i]
        copynumber_DATA <- DATA_cancer_type_cn[[i]]
        #   Find statistics for each cancer type
        DATA_statistics <- get_WGD_stats_from_data(
            copynumber_DATA, DATA_wgd, copynumber_coordinates, cn_info,
            list_targets = c("WGD_proportion", "FGA_difference", "WGD_FGA_by_sample")
        )
        #   Get WGD proportion for each cancer type
        WGD_proportion <- DATA_statistics$WGD_proportion
        #   Get WGD-associated FGA difference for each cancer type
        WGD_increased_FGA <- DATA_statistics$FGA_difference
        #   Get table of sample WGD and FGA for each cancer type
        WGD_FGA_by_sample <- DATA_statistics$WGD_FGA_by_sample
        df_tmp <- data.frame(
            cancer_type = cancer_type,
            WGD = WGD_FGA_by_sample$WGD,
            FGA = WGD_FGA_by_sample$FGA
        )
        if (i == 1) {
            df_WGD_FGA_by_sample <- df_tmp
        } else {
            df_WGD_FGA_by_sample <- rbind(df_WGD_FGA_by_sample, df_tmp)
        }
        df_WGD_FGA[i, ] <- c(cancer_type, WGD_proportion, WGD_increased_FGA)
    }
    cat("\n")
    df_WGD_FGA$WGD <- as.numeric(df_WGD_FGA$WGD)
    df_WGD_FGA$WGD_increased_FGA <- as.numeric(df_WGD_FGA$WGD_increased_FGA)
    write.csv(df_WGD_FGA, file = paste0(plotname, ".csv"))
    if (length(which(df_WGD_FGA$WGD == 0)) > 0) {
        df_WGD_FGA <- df_WGD_FGA[-which(df_WGD_FGA$WGD == 0), ]
    }
    #---------------------------Find WGD proportions in each cancer type
    DATA_wgd_proportion <- rep(0, length(DATA_cancer_types))
    for (i in 1:length(DATA_cancer_types)) {
        mini_DATA_wgd <- DATA_wgd[which(DATA_wgd$samplename %in% DATA_cancer_type_sample_ids[[i]] & DATA_wgd$wgd_uncertain == FALSE), ]
        DATA_wgd_proportion[i] <- 100 * length(which(mini_DATA_wgd$wgd_status == "wgd")) / length(mini_DATA_wgd$wgd_status)
    }
    #------------Find count and strength of TSG arms in each cancer type
    FIT_tsg_count <- rep(0, length(DATA_cancer_types))
    FIT_tsg_mean_selection_rate <- rep(0, length(DATA_cancer_types))
    FIT_tsg_max_selection_rate <- rep(0, length(DATA_cancer_types))
    FIT_onc_count <- rep(0, length(DATA_cancer_types))
    FIT_onc_mean_selection_rate <- rep(0, length(DATA_cancer_types))
    FIT_onc_max_selection_rate <- rep(0, length(DATA_cancer_types))
    for (i in 1:length(DATA_cancer_types)) {
        cancer_types_fit <- read.csv(paste0(DATA_cancer_types[i], "_fitted_parameters.csv"), header = TRUE)
        cancer_types_fit_selection_rates <- cancer_types_fit$Best_value[which(cancer_types_fit$Type == "Arm_selection_rate")]
        FIT_tsg_count[i] <- length(which(cancer_types_fit_selection_rates < 1))
        if (FIT_tsg_count[i] == 0) {
            FIT_tsg_mean_selection_rate[i] <- 1
            FIT_tsg_max_selection_rate[i] <- 1
        } else {
            tmp <- cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates < 1)]
            FIT_tsg_mean_selection_rate[i] <- 1 / prod(tmp)^(1 / length(tmp))
            FIT_tsg_max_selection_rate[i] <- 1 / max(tmp)
            # FIT_tsg_mean_selection_rate[i] <- 1 / mean(cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates < 1)])
            # FIT_tsg_max_selection_rate[i] <- 1 / max(cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates < 1)])
        }
        FIT_onc_count[i] <- length(which(cancer_types_fit_selection_rates > 1))
        if (FIT_onc_count[i] == 0) {
            FIT_onc_mean_selection_rate[i] <- 1
            FIT_onc_max_selection_rate[i] <- 1
        } else {
            tmp <- cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates > 1)]
            FIT_onc_mean_selection_rate[i] <- prod(tmp)^(1 / length(tmp))
            FIT_onc_max_selection_rate[i] <- max(tmp)
            # FIT_onc_mean_selection_rate[i] <- mean(cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates > 1)])
            # FIT_onc_max_selection_rate[i] <- max(cancer_types_fit_selection_rates[which(cancer_types_fit_selection_rates > 1)])
        }
    }

    #-------------Plot data FGA violin plots by cancer type & WGD status
    geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                                  draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                                  show.legend = NA, inherit.aes = TRUE) {
        layer(
            data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
            position = position, show.legend = show.legend, inherit.aes = inherit.aes,
            params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
        )
    }


    #   Function to make splitted violin plots
    GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
        draw_group = function(self, data, ..., draw_quantiles = NULL) {
            data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
            grp <- data[1, "group"]
            newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
            newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
            newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

            if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                    1))
                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                aesthetics$alpha <- rep(1, nrow(quantiles))
                both <- cbind(quantiles, aesthetics)
                quantile_grob <- GeomPath$draw_panel(both, ...)
                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
            } else {
                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
            }
        }
    )
    filename <- paste0(plotname, "FGA_from_data.jpeg")
    jpeg(filename, width = 2000, height = 1100)
    p <- ggplot(df_WGD_FGA_by_sample, aes(x = cancer_type, y = FGA, fill = WGD)) +
        geom_violin(scale = "width")
    # geom_split_violin(width = 1, alpha = 0.2, scale = "width")
    print(p)
    dev.off()
    print(df_WGD_FGA_by_sample)

    #---------------------------Plot relationship between WGD proportion
    #-----------------------------------------------and aneuploidy score
    filename <- paste0(plotname, "_WGD_vs_FGA_DIFFERENCE_from_data.jpeg")
    jpeg(filename, width = 2000, height = 1100)
    p <- ggplot(df_WGD_FGA, aes(x = WGD, y = WGD_increased_FGA)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_type), size = 10, box.padding = 1) +
        # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        # ylim(0, NA) +
        xlab("WGD proportion") +
        ylab("FGA difference") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #-------------------Plot relationship between fitted selection rates
    #---------------------------------------and observed WGD proportions
    #---Prepare the dataframe for plotting
    df_plot <- data.frame(
        cancer_types = DATA_cancer_types,
        WGD = DATA_wgd_proportion,
        TSG_count = FIT_tsg_count,
        TSG_mean_selection_rate = FIT_tsg_mean_selection_rate,
        TSG_max_selection_rate = FIT_tsg_max_selection_rate,
        ONC_count = FIT_onc_count,
        ONC_mean_selection_rate = FIT_onc_mean_selection_rate,
        ONC_max_selection_rate = FIT_onc_max_selection_rate
    )
    #---Find p-values for correlation between WGD status and either count or mean selection rate of TSG/ONC arms
    tmp <- cor.test(FIT_tsg_count, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    p_val_TSG_count <- tmp$p.value
    tmp <- cor.test(FIT_tsg_mean_selection_rate, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    p_val_TSG_mean_selection_rate <- tmp$p.value
    tmp <- cor.test(FIT_onc_count, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    p_val_ONC_count <- tmp$p.value
    tmp <- cor.test(FIT_onc_mean_selection_rate, DATA_wgd_proportion, method = "spearman", exact = FALSE)
    p_val_ONC_mean_selection_rate <- tmp$p.value
    #---Positions for p-values
    x_left <- 0.0
    x_right <- 0.75
    y_down <- 0.0
    y_up <- 0.95
    #---Plot relationship between WGD status and count of GAIN/LOSS arms
    filename <- paste0(plotname, "_WGD_vs_GAIN&LOSS_counts.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    p <- ggplot(df_plot, aes(x = TSG_count, y = ONC_count, color = WGD)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_types), size = 10, box.padding = 1) +
        annotate("segment",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            xend = min(df_plot$TSG_count) + (x_left + 0.05) * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$ONC_count) + y_up * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            yend = min(df_plot$ONC_count) + y_up * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_count) + (x_left + 0.07) * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$ONC_count) + y_up * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            label = paste0("p.val=", scientific(p_val_TSG_count, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        annotate("segment",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            xend = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$ONC_count) + y_up * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            yend = min(df_plot$ONC_count) + (y_up + 0.05) * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$ONC_count) + (y_up + 0.07) * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            label = paste0("p.val=", scientific(p_val_ONC_count, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        xlab("Count of LOSS arms") +
        ylab("Count of GAIN arms") +
        labs(fill = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #---Plot relationship between WGD status and selection rates of GAIN/LOSS arms
    filename <- paste0(plotname, "_WGD_vs_GAIN&LOSS_selection_rates.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    p <- ggplot(df_plot, aes(x = TSG_mean_selection_rate, y = ONC_mean_selection_rate, color = WGD)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_types), size = 10, box.padding = 1) +
        annotate("segment",
            x = min(df_plot$TSG_mean_selection_rate) + x_left * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            xend = min(df_plot$TSG_mean_selection_rate) + (x_left + 0.05) * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            yend = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_mean_selection_rate) + (x_left + 0.07) * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_TSG_mean_selection_rate, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        annotate("segment",
            x = min(df_plot$TSG_mean_selection_rate) + x_left * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            xend = min(df_plot$TSG_mean_selection_rate) + x_left * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            yend = min(df_plot$ONC_mean_selection_rate) + (y_up + 0.05) * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_mean_selection_rate) + x_left * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            y = min(df_plot$ONC_mean_selection_rate) + (y_up + 0.07) * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_ONC_mean_selection_rate, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        xlab("Mean selection rate of LOSS arms") +
        ylab("Mean selection rate of GAIN arms") +
        labs(fill = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #---Plot relationship between WGD status and count/mean selection rate of LOSS arms
    filename <- paste0(plotname, "_WGD_vs_LOSS_fitted.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    p <- ggplot(df_plot, aes(x = TSG_count, y = TSG_mean_selection_rate, color = WGD)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_types), size = 10, box.padding = 1) +
        annotate("segment",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            xend = min(df_plot$TSG_count) + (x_left + 0.05) * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$TSG_mean_selection_rate) + y_up * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            yend = min(df_plot$TSG_mean_selection_rate) + y_up * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_count) + (x_left + 0.07) * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$TSG_mean_selection_rate) + y_up * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_TSG_count, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        annotate("segment",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            xend = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$TSG_mean_selection_rate) + y_up * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            yend = min(df_plot$TSG_mean_selection_rate) + (y_up + 0.05) * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$TSG_count) + x_left * (max(df_plot$TSG_count) - min(df_plot$TSG_count)),
            y = min(df_plot$TSG_mean_selection_rate) + (y_up + 0.07) * (max(df_plot$TSG_mean_selection_rate) - min(df_plot$TSG_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_TSG_mean_selection_rate, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        xlab("Count of LOSS arms") +
        ylab("Mean selection rate of LOSS arms") +
        labs(fill = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #---Plot relationship between WGD status and count/mean selection rate of GAIN arms
    filename <- paste0(plotname, "_WGD_vs_GAIN_fitted.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    p <- ggplot(df_plot, aes(x = ONC_count, y = ONC_mean_selection_rate, color = WGD)) +
        geom_point(size = 10) +
        geom_text_repel(aes(label = cancer_types), size = 10, box.padding = 1, point.padding = 0.5) +
        annotate("segment",
            x = min(df_plot$ONC_count) + x_left * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            xend = min(df_plot$ONC_count) + (x_left + 0.05) * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            yend = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$ONC_count) + (x_left + 0.07) * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_ONC_count, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        annotate("segment",
            x = min(df_plot$ONC_count) + x_left * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            xend = min(df_plot$ONC_count) + x_left * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            y = min(df_plot$ONC_mean_selection_rate) + y_up * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            yend = min(df_plot$ONC_mean_selection_rate) + (y_up + 0.05) * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            colour = "black", size = 2, alpha = 1, arrow = arrow()
        ) +
        annotate("text",
            x = min(df_plot$ONC_count) + x_left * (max(df_plot$ONC_count) - min(df_plot$ONC_count)),
            y = min(df_plot$ONC_mean_selection_rate) + (y_up + 0.07) * (max(df_plot$ONC_mean_selection_rate) - min(df_plot$ONC_mean_selection_rate)),
            label = paste0("p.val=", scientific(p_val_ONC_mean_selection_rate, digits = 3)), colour = "black", size = 8, hjust = 0
        ) +
        xlab("Count of GAIN arms") +
        ylab("Mean selection rate of GAIN arms") +
        labs(fill = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    return(df_WGD_FGA)
}
