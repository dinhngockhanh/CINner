#----------------------Function to assign parameters to proper positions
bulk_gene_assign_paras <- function(model_variables, parameter_IDs, parameters) {
    for (i in 1:length(parameter_IDs)) {
        parameter_ID <- parameter_IDs[i]
        if (parameter_ID %in% model_variables$general_variables$Variable) {
            model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameters[i]
        } else if (parameter_ID %in% model_variables$driver_library$Gene_ID) {
            model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == parameter_ID)] <- parameters[i]
        }
    }
    model_variables <- BUILD_driver_library(
        model_variables = model_variables,
        table_gene_selection_rates = model_variables$driver_library
    )
    return(model_variables)
}
#------Function to extract simulation gene-level mut/amp/del frequencies
bulk_gene_get_gene_counts <- function(SIMS_chromosome, list_targets, driver_library) {
    #-----------------------Get mut/amp/del frequencies from simulations
    count_mut <- matrix(0, nrow = length(SIMS_chromosome), ncol = length(list_targets))
    count_gain <- matrix(0, nrow = length(SIMS_chromosome), ncol = length(list_targets))
    count_loss <- matrix(0, nrow = length(SIMS_chromosome), ncol = length(list_targets))
    for (iteration in 1:length(SIMS_chromosome)) {
        simulation <- SIMS_chromosome[[iteration]]
        all_sample_genotype <- simulation$sample$all_sample_genotype
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
        max_genotype_profile <- sample_genotype_unique_profile[[max_genotype]]
        max_genotype_drivers <- sample_genotype_unique_drivers[[max_genotype]]
        #   Get the genotype's driver mutations
        if (nrow(max_genotype_drivers) > 0) {
            for (i in 1:nrow(max_genotype_drivers)) {
                loc <- which(list_targets == driver_library$Gene_ID[max_genotype_drivers[i, 1]])
                count_mut[iteration, loc] <- 1
            }
        }
        #   Get the genotype's driver amplifications and deletions
        for (i in 1:length(list_targets)) {
            loc <- which(driver_library$Gene_ID == list_targets[i])
            chrom <- driver_library$Chromosome[loc]
            bin <- driver_library$Bin[loc]
            cn <- max_genotype_profile$copy[which(max_genotype_profile$chr == chrom)][bin]
            if (cn > 2) count_gain[iteration, i] <- 1
            if (cn < 2) count_loss[iteration, i] <- 1
        }
    }
    #------------------------------------------Output frequency matrices
    output <- list()
    output$count_mut <- count_mut
    output$count_gain <- count_gain
    output$count_loss <- count_loss
    return(output)
}
#---------Function to extract average gene-level mut/amp/del frequencies
bulk_gene_get_gene_freqs <- function(SIMS_chromosome, list_targets, driver_library) {
    #-----------------------Get mut/amp/del frequencies from simulations
    count_matrices <- bulk_gene_get_gene_counts(SIMS_chromosome, list_targets, driver_library)
    count_mut <- count_matrices$count_mut
    count_gain <- count_matrices$count_gain
    count_loss <- count_matrices$count_loss
    #-----------Get average mut/amp/del frequencies from all simulations
    gene_freqs <- c(
        colMeans(count_mut),
        colMeans(count_gain),
        colMeans(count_loss)
    )
    return(gene_freqs)
}
#-------------------------------------Objective function for ABC fitting
bulk_gene_func_ABC <- function(parameters, parameter_IDs, model_variables, list_targets) {
    #   Assign parameters in model variables
    model_variables <- bulk_gene_assign_paras(model_variables, parameter_IDs, parameters)
    driver_library <- model_variables$driver_library
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE, build_cn = FALSE,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile", "sample_genotype_unique_drivers"),
        R_libPaths = R_libPaths
    )
    #   Statistics = gene-level mutation/amplification/deletion delta
    SIMS_gene_freqs <- bulk_gene_get_gene_freqs(SIMS_chromosome, list_targets, driver_library)
    stat <- SIMS_gene_freqs
}
#-----------------Function to choose one best parameter from a posterior
bulk_gene_get_best_para <- function(data_rf, model_rf, obs_rf, post_rf) {
    df_dist <- densityPlot_df(model_rf, obs_rf, data_rf)
    best_para <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
}

#' @export
library_bulk_gene <- function(library_name,
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
    bulk_gene_func_ABC <<- bulk_gene_func_ABC
    bulk_gene_assign_paras <<- bulk_gene_assign_paras
    bulk_gene_get_gene_freqs <<- bulk_gene_get_gene_freqs
    bulk_gene_get_gene_counts <<- bulk_gene_get_gene_counts
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c(
        "n_samples", "list_targets", "sim_param", "parameter_IDs", "model_variables", "gainloss_SIMS",
        "bulk_gene_func_ABC", "bulk_gene_assign_paras", "bulk_gene_get_gene_freqs", "bulk_gene_get_gene_counts", "R_libPaths",
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
        stat <- bulk_gene_func_ABC(parameters, parameter_IDs, model_variables, list_targets)
        return(stat)
    })
    stopCluster(cl)
    #   Group simulated statistics into one table
    sim_stat <- matrix(0, nrow = ABC_simcount, ncol = 3 * length(list_targets))
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
fitting_bulk_gene <- function(library_name,
                              model_name,
                              sample_ids_DATA,
                              driver_events_DATA,
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
    library(dplyr)
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
    #-------------Find mut/amp/del frequencies for entire driver library
    #   Find driver event count matrices
    DATA_count_mut <- matrix(0, nrow = length(sample_ids_DATA), ncol = length(list_targets))
    DATA_count_gain <- matrix(0, nrow = length(sample_ids_DATA), ncol = length(list_targets))
    DATA_count_loss <- matrix(0, nrow = length(sample_ids_DATA), ncol = length(list_targets))
    for (i in 1:length(sample_ids_DATA)) {
        sample_id <- sample_ids_DATA[i]
        sample_driver_events_DATA <- driver_events_DATA[which(driver_events_DATA$sample_id == sample_id), ]
        if (nrow(sample_driver_events_DATA) == 0) next
        for (j in 1:nrow(sample_driver_events_DATA)) {
            gene <- sample_driver_events_DATA$gene[j]
            type <- sample_driver_events_DATA$top_category[j]
            category <- sample_driver_events_DATA$category[j]
            if (!gene %in% list_targets) next
            if (type == "mutational") DATA_count_mut[i, which(list_targets == gene)] <- 1
            if (type == "CNA") {
                if (category == "coding_amplification") {
                    DATA_count_gain[i, which(list_targets == gene)] <- 1
                } else if (category == "coding_deletion") {
                    DATA_count_loss[i, which(list_targets == gene)] <- 1
                }
            }
        }
    }
    #   Find statistics target
    DATA_target <- c(
        colMeans(DATA_count_mut),
        colMeans(DATA_count_gain),
        colMeans(DATA_count_loss)
    )
    # ====================================FITTING WITH ABC RANDOM FOREST
    #--------------------------------------------Fit parameters with ABC
    #---Dataframe for data observation
    obs_rf <- data.frame(matrix(DATA_target, nrow = 1))
    colnames(obs_rf) <- paste0("freqs_", 1:ncol(obs_rf))
    #---Dataframe for parameters for reference
    all_paras <- data.frame(sim_param)
    colnames(all_paras) <- parameter_IDs
    #---Dataframe for corresponding statistics for reference
    all_data <- data.frame(sim_stat)
    colnames(all_data) <- paste0("freqs_", 1:ncol(all_data))
    #---Fit each parameter with ABC-rf
    layout <- matrix(NA, nrow = 7, ncol = ceiling(length(parameter_IDs) / 7))
    gs <- list()
    id <- 0
    for (para in 1:nrow(list_parameters)) {
        para_ID <- list_parameters$Variable[para]
        cat(paste("ABC for parameter ", para_ID, " [", para, "/", nrow(list_parameters), "]", "\n", sep = ""))
        #   Train the random forest
        data_rf <- cbind(all_paras[para_ID], all_data)
        colnames(data_rf)[1] <- "para"
        f <- as.formula("para ~.")
        model_rf <- regAbcrf(formula = f, data_rf, paral = TRUE, ncores = n_cores)
        #   Predict posterior distribution based on found random forest
        post_rf <- predict(model_rf, obs_rf, data_rf, paral = TRUE, ncores = n_cores)
        #   Choose best value from posterior distribution
        best_rf <- bulk_gene_get_best_para(data_rf, model_rf, obs_rf, post_rf)
        #   Save results for fitting this parameter
        ABC_output <- list()
        ABC_output$para_ID <- para_ID
        ABC_output$data_rf <- data_rf
        ABC_output$model_rf <- model_rf
        ABC_output$obs_rf <- obs_rf
        ABC_output$post_rf <- post_rf
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
            model_rf, obs_rf, data_rf,
            protocol = "arm",
            fontsize = 20,
            chosen_para = best_rf,
            color_prior = "lightblue", color_posterior = "darkblue", color_vline = "blue",
            main = para_ID
        )
    }
    #   Plot the prior, posterior and chosen best parameter for all variables
    filename <- paste0(model_name, "_ABC_all.jpeg")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
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
    #-----------------Analysis of fitted driver frequencies against data
    #   Assign parameters in model variables
    model_variables <- bulk_gene_assign_paras(model_variables, list_parameters$Variable, list_parameters$Best_value)
    #   Make simulations using best parameters
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_samples, stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE, build_cn = FALSE, compute_parallel = TRUE, n_cores = n_cores,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile", "sample_genotype_unique_drivers"),
        R_libPaths = R_libPaths
    )
    #   Find driver event count matrices
    driver_library <- model_variables$driver_library
    SIMS_best <- bulk_gene_get_gene_counts(SIMS_chromosome, list_targets, driver_library)
    SIMS_count_mut <- SIMS_best$count_mut
    SIMS_count_gain <- SIMS_best$count_gain
    SIMS_count_loss <- SIMS_best$count_loss
    #   Compute averages of frequencies
    DATA_mean_freqs <- c(100 * colMeans(DATA_count_mut), 100 * colMeans(DATA_count_gain), 100 * colMeans(DATA_count_loss))
    SIMS_mean_freqs <- c(100 * colMeans(SIMS_count_mut), 100 * colMeans(SIMS_count_gain), 100 * colMeans(SIMS_count_loss))
    #   Compute standard errors of frequencies from simulations
    N_bootstrap <- 10000
    DATA_stat_bootstrap <- matrix(0, nrow = N_bootstrap, ncol = length(DATA_mean_freqs))
    SIMS_stat_bootstrap <- matrix(0, nrow = N_bootstrap, ncol = length(SIMS_mean_freqs))
    for (i in 1:N_bootstrap) {
        DATA_locs <- sample(1:length(sample_ids_DATA), length(sample_ids_DATA), replace = TRUE)
        DATA_count_mut_tmp <- DATA_count_mut[DATA_locs, ]
        DATA_count_gain_tmp <- DATA_count_gain[DATA_locs, ]
        DATA_count_loss_tmp <- DATA_count_loss[DATA_locs, ]
        DATA_stat_bootstrap[i, ] <- c(100 * colMeans(DATA_count_mut_tmp), 100 * colMeans(DATA_count_gain_tmp), 100 * colMeans(DATA_count_loss_tmp))

        SIMS_locs <- sample(1:n_samples, n_samples, replace = TRUE)
        SIMS_count_mut_tmp <- SIMS_count_mut[SIMS_locs, ]
        SIMS_count_gain_tmp <- SIMS_count_gain[SIMS_locs, ]
        SIMS_count_loss_tmp <- SIMS_count_loss[SIMS_locs, ]
        SIMS_stat_bootstrap[i, ] <- c(100 * colMeans(SIMS_count_mut_tmp), 100 * colMeans(SIMS_count_gain_tmp), 100 * colMeans(SIMS_count_loss_tmp))
    }
    DATA_sd_freqs <- apply(DATA_stat_bootstrap, 2, sd)
    SIMS_sd_freqs <- apply(SIMS_stat_bootstrap, 2, sd)
    #   Make dataframe of data frequencies and simulated frequencies
    Gene_ID <- list_targets
    df_freq <- data.frame(
        Gene_ID = rep(list_targets, 6),
        Freq = c(DATA_mean_freqs, SIMS_mean_freqs),
        Type = rep(c(rep("Mutation", length(list_targets)), rep("Gain", length(list_targets)), rep("Loss", length(list_targets))), 2),
        Source = c(rep("Data", length(DATA_mean_freqs)), rep("Simulation", length(SIMS_mean_freqs))),
        SD = c(DATA_sd_freqs, SIMS_sd_freqs)
    )
    #   Plot comparison vs data
    filename <- paste0(model_name, "_against_data.jpeg")
    jpeg(file = filename, width = 2000, height = 1100)
    p <- ggplot(df_freq, aes(x = Gene_ID, y = Freq, fill = Type, color = Type, Source = Source)) +
        geom_bar(stat = "identity", width = 0.5, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Simulation"), alpha = 0.6) +
        geom_errorbar(aes(ymin = Freq - SD, ymax = Freq + SD), width = 0.5, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Simulation")) +
        geom_point(size = 5, shape = 21, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Data")) +
        geom_errorbar(aes(ymin = Freq - SD, ymax = Freq + SD), width = 0.5, position = position_dodge(width = 0.5), data = df_freq %>% filter(Source == "Data")) +
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
    #------------Analysis of fitted parameters against mut/amp/del freqs
    #   Prepare plot dataframe with fitted selection rates and mut/amp/del freqs
    plot_table <- data.frame(Gene = list_targets)
    plot_table$Mut_freq_spec <- 0
    plot_table$Gain_freq_spec <- 0
    plot_table$Loss_freq_spec <- 0
    plot_table$Selection_rate <- 0
    for (row in 1:nrow(plot_table)) {
        gene <- plot_table$Gene[row]
        plot_table$Mut_freq_spec[row] <- 100 * DATA_target[which(list_targets == gene)]
        plot_table$Gain_freq_spec[row] <- 100 * DATA_target[which(list_targets == gene) + length(list_targets)]
        plot_table$Loss_freq_spec[row] <- 100 * DATA_target[which(list_targets == gene) + 2 * length(list_targets)]
        plot_table$Selection_rate[row] <- parameters_best[which(list_parameters$Variable == gene)]
    }
    #   Plot selection rates vs cancer-specific mut frequencies
    filename <- paste(model_name, "_selection_rate_analysis.jpeg", sep = "")
    jpeg(filename, width = 1000, height = 1100)
    tmp <- cor.test(plot_table$Selection_rate, plot_table$Mut_freq_spec)
    Pearson_r <- tmp$estimate
    Pearson_p_val <- tmp$p.value
    grob1 <- grobTree(
        textGrob(
            paste("r=", round(Pearson_r, 3), ", p-value=", formatC(Pearson_p_val, format = "e", digits = 2)),
            x = 0.1, y = 0.9, just = "left",
            gp = gpar(col = "blue", fontsize = 30, fontface = "bold"),
        )
    )
    p <- ggplot(plot_table, aes(x = Selection_rate, y = Mut_freq_spec)) +
        geom_point(size = 30, color = "blue") +
        geom_smooth(method = lm, color = "blue") +
        annotation_custom(grob1) +
        geom_text(aes(label = Gene), size = 10, color = "white") +
        xlab("Selection rates") +
        ylab("Mutation frequencies") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 40)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    print(p)
    dev.off()
}
