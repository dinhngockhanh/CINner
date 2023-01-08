#' @export
simulator_multivar <- function(model_variables_base = list(),
                               model_prefix = "",
                               var1_name = "",
                               var1_vals = c(),
                               var2_name = "",
                               var2_vals = c(),
                               extra_var = NULL,
                               n_simulations = 0,
                               stage_final = 3,
                               n_clones_min = 0,
                               n_clones_max = Inf,
                               save_simulation = TRUE,
                               neutral_variations = FALSE,
                               internal_nodes_cn_info = FALSE,
                               save_newick_tree = FALSE,
                               save_cn_profile = FALSE,
                               save_cn_clones = FALSE,
                               format_cn_profile = "long",
                               model_readcount = FALSE,
                               model_readcount_base = "all",
                               pseudo_corrected_readcount = FALSE,
                               HMM = FALSE,
                               HMM_containner = "docker",
                               folder_workplace = NULL,
                               report_progress = TRUE,
                               compute_parallel = TRUE,
                               seed = Inf,
                               output_variables = c(),
                               n_cores = NULL,
                               R_libPaths = NULL,
                               plot = FALSE) {
    library(scales)
    if (is.vector(var1_vals)) {
        rows <- length(var1_vals):1
    } else if (is.matrix(var1_vals)) {
        rows <- ncol(var1_vals):1
    }
    if (is.vector(var2_vals)) {
        cols <- length(var2_vals):1
    } else if (is.matrix(var2_vals)) {
        cols <- ncol(var2_vals):1
    }
    ind <- 0
    for (row in rows) {
        for (col in cols) {
            ind <- ind + 1
            #-----------------------Create model variables for the batch
            #   Initialize parameter set
            model_variables <- model_variables_base
            #   Fix variable 1 in parameter set
            if (var1_name == "prob_CN_whole_genome_duplication") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var1_vals[row]
            } else if (var1_name == "prob_CN_missegregation") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- var1_vals[row]
            } else if (var1_name == "delta_selection") {
                delta_sel_gain_per_arm <- sqrt(1 + var1_vals[1, row]) - 1
                delta_sel_loss_per_arm <- sqrt(1 + var1_vals[2, row]) - 1
                chrom_gains <- sample(model_variables_base$cn_info$Chromosome, round(length(model_variables_base$cn_info$Chromosome) / 2), replace = FALSE)
                chrom_losses <- setdiff(model_variables_base$cn_info$Chromosome, chrom_gains)
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_gains)] <- 1 + delta_sel_gain_per_arm
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_losses)] <- 1 / (1 + delta_sel_loss_per_arm)
            } else if (var1_name == "bound_homozygosity") {
                model_variables$selection_model$Value[which(model_variables$selection_model$Variable == "bound_homozygosity")] <- var1_vals[row]
            }
            #   Fix variable 2 in parameter set
            if (var2_name == "prob_CN_whole_genome_duplication") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var2_vals[col]
            } else if (var2_name == "prob_CN_missegregation") {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- var2_vals[col]
            } else if (var2_name == "delta_selection") {
                delta_sel_gain_per_arm <- sqrt(1 + var2_vals[1, col]) - 1
                delta_sel_loss_per_arm <- sqrt(1 + var2_vals[2, col]) - 1
                chrom_gains <- sample(model_variables_base$cn_info$Chromosome, round(length(model_variables_base$cn_info$Chromosome) / 2), replace = FALSE)
                chrom_losses <- setdiff(model_variables_base$cn_info$Chromosome, chrom_gains)
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_gains)] <- 1 + delta_sel_gain_per_arm
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_losses)] <- 1 / (1 + delta_sel_loss_per_arm)
            } else if (var2_name == "bound_homozygosity") {
                model_variables$selection_model$Value[which(model_variables$selection_model$Variable == "bound_homozygosity")] <- var2_vals[col]
            }





            #   Fix additional variables if necessary
            if (!is.null(extra_var)) {
                for (i in 1:ceiling(length(extra_var) / 2)) {
                    variable_type <- extra_var[2 * i - 1]
                    if (variable_type == "prob_neutral_CN_missegregation/prob_CN_missegregation") {
                        model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_neutral_CN_missegregation")] <- as.numeric(extra_var[2 * i]) * as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")])
                    }
                }
            }





            #   Save model variables
            if (var1_name == "delta_selection") {
                model_name <- paste(model_prefix, "_", var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_", var2_name, "=", scientific(var2_vals[col]), sep = "")
            } else if (var2_name == "delta_selection") {
                model_name <- paste(model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]), sep = "")
            } else {
                model_name <- paste(model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[col]), sep = "")
            }
            SAVE_model_variables(model_name = model_name, model_variables = model_variables)
            #   Create simulations
            cat("=======================================================\n")
            cat("=======================================================\n")
            cat("=======================================================\n")
            cat(paste("\nSIMULATIONS FOR BATCH ", ind, "/", length(rows) * length(cols), "...\n", sep = ""))
            if (var1_name == "delta_selection") {
                cat(paste(var1_name, " = ", scientific(var1_vals[1, row]), " & ", scientific(var1_vals[2, row]), "\n", sep = ""))
            } else {
                cat(paste(var1_name, " = ", scientific(var1_vals[row]), "\n", sep = ""))
            }
            if (var2_name == "delta_selection") {
                cat(paste(var2_name, " = ", scientific(var2_vals[1, col]), " & ", scientific(var2_vals[2, col]), "\n", sep = ""))
            } else {
                cat(paste(var2_name, " = ", scientific(var2_vals[col]), "\n", sep = ""))
            }
            #---------------------------Create simulations for the batch
            start_time <- Sys.time()
            tmp <- simulator_full_program(
                model = model_name,
                n_simulations = n_simulations,
                stage_final = stage_final,
                n_clones_min = n_clones_min,
                n_clones_max = n_clones_max,
                neutral_variations = neutral_variations,
                internal_nodes_cn_info = internal_nodes_cn_info,
                save_newick_tree = save_newick_tree,
                save_cn_profile = save_cn_profile,
                save_cn_clones = save_cn_clones,
                format_cn_profile = format_cn_profile,
                model_readcount = model_readcount,
                model_readcount_base = model_readcount_base,
                pseudo_corrected_readcount = pseudo_corrected_readcount,
                HMM = HMM,
                HMM_containner = HMM_containner,
                folder_workplace = folder_workplace,
                report_progress = report_progress,
                compute_parallel = compute_parallel,
                seed = seed,
                output_variables = output_variables,
                n_cores = n_cores,
                save_simulation = save_simulation,
                R_libPaths = R_libPaths
            )
            end_time <- Sys.time()
            print(end_time - start_time)
            cat("\n")
            #-----------------------------Plot simulations for the batch
            if (plot == TRUE) {
                plot_all(
                    model = model_name,
                    n_simulations = n_simulations,
                    unit_time = model_variables$general_variables$Unit[which(model_variables$general_variables$Variable == "T_end_time")],
                    folder_workplace = folder_workplace,
                    folder_plots = folder_workplace,
                    compute_parallel = compute_parallel,
                    R_libPaths = R_libPaths
                )
            }
            #-----Move model variable files for the batch into workplace
            if (file.exists(paste0(model_name, "-input-copy-number-blocks.csv"))) {
                file.copy(from = paste0(model_name, "-input-copy-number-blocks.csv"), to = paste0(folder_workplace, "/", model_name, "-input-copy-number-blocks.csv"))
                file.remove(paste0(model_name, "-input-copy-number-blocks.csv"))
            }
            if (file.exists(paste0(model_name, "-input-gc.csv"))) {
                file.copy(from = paste0(model_name, "-input-gc.csv"), to = paste0(folder_workplace, "/", model_name, "-input-gc.csv"))
                file.remove(paste0(model_name, "-input-gc.csv"))
            }
            if (file.exists(paste0(model_name, "-input-initial-cn-profiles.csv"))) {
                file.copy(from = paste0(model_name, "-input-initial-cn-profiles.csv"), to = paste0(folder_workplace, "/", model_name, "-input-initial-cn-profiles.csv"))
                file.remove(paste0(model_name, "-input-initial-cn-profiles.csv"))
            }
            if (file.exists(paste0(model_name, "-input-initial-others.csv"))) {
                file.copy(from = paste0(model_name, "-input-initial-others.csv"), to = paste0(folder_workplace, "/", model_name, "-input-initial-others.csv"))
                file.remove(paste0(model_name, "-input-initial-others.csv"))
            }
            if (file.exists(paste0(model_name, "-input-population-dynamics.csv"))) {
                file.copy(from = paste0(model_name, "-input-population-dynamics.csv"), to = paste0(folder_workplace, "/", model_name, "-input-population-dynamics.csv"))
                file.remove(paste0(model_name, "-input-population-dynamics.csv"))
            }
            if (file.exists(paste0(model_name, "-input-sampling.csv"))) {
                file.copy(from = paste0(model_name, "-input-sampling.csv"), to = paste0(folder_workplace, "/", model_name, "-input-sampling.csv"))
                file.remove(paste0(model_name, "-input-sampling.csv"))
            }
            if (file.exists(paste0(model_name, "-input-selection-chrom-arm.csv"))) {
                file.copy(from = paste0(model_name, "-input-selection-chrom-arm.csv"), to = paste0(folder_workplace, "/", model_name, "-input-selection-chrom-arm.csv"))
                file.remove(paste0(model_name, "-input-selection-chrom-arm.csv"))
            }
            if (file.exists(paste0(model_name, "-input-selection-genes.csv"))) {
                file.copy(from = paste0(model_name, "-input-selection-genes.csv"), to = paste0(folder_workplace, "/", model_name, "-input-selection-genes.csv"))
                file.remove(paste0(model_name, "-input-selection-genes.csv"))
            }
            if (file.exists(paste0(model_name, "-input-selection-model.csv"))) {
                file.copy(from = paste0(model_name, "-input-selection-model.csv"), to = paste0(folder_workplace, "/", model_name, "-input-selection-model.csv"))
                file.remove(paste0(model_name, "-input-selection-model.csv"))
            }
            if (file.exists(paste0(model_name, "-input-variables.csv"))) {
                file.copy(from = paste0(model_name, "-input-variables.csv"), to = paste0(folder_workplace, "/", model_name, "-input-variables.csv"))
                file.remove(paste0(model_name, "-input-variables.csv"))
            }
        }
    }
}

#' @export
statistics_multivar <- function(model_prefix = "",
                                folder_workplace = "",
                                var1_name = "",
                                var1_vals = c(),
                                var2_name = "",
                                var2_vals = c(),
                                n_simulations = 0,
                                compute_parallel = TRUE,
                                n_cores = NULL,
                                R_libPaths = NULL) {
    library(data.table)
    library(scales)
    library(ggplot2)
    library(viridis)
    library(hrbrthemes)
    library(scatterpie)
    if (var1_name == "prob_CN_whole_genome_duplication") {
        var1_lab <- "Probability of WGD"
    } else if (var1_name == "delta_selection") {
        var1_lab <- "Selection rate"
    } else if (var1_name == "prob_CN_missegregation") {
        var1_lab <- "Probability of missegregation"
    }
    if (var2_name == "prob_CN_whole_genome_duplication") {
        var2_lab <- "Probability of WGD"
    } else if (var2_name == "delta_selection") {
        var2_lab <- "Selection rate"
    } else if (var2_name == "prob_CN_missegregation") {
        var2_lab <- "Probability of missegregation"
    }
    if (is.vector(var1_vals)) {
        rows <- length(var1_vals):1
    } else if (is.matrix(var1_vals)) {
        rows <- ncol(var1_vals):1
    }
    if (is.vector(var2_vals)) {
        cols <- length(var2_vals):1
    } else if (is.matrix(var2_vals)) {
        cols <- ncol(var2_vals):1
    }
    #-----------------------------Whether to plot WGD-related statistics
    if ((var1_name == "prob_CN_whole_genome_duplication") | (var2_name == "prob_CN_whole_genome_duplication")) {
        plot_WGD <- TRUE
    } else {
        plot_WGD <- FALSE
    }
    #------------Create dataframe for all statistics for all simulations
    df_stats_simulations <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_stats_simulations) <- c("var1", "var2", "sim", "stat", "val")







    # #------------------------------------Get statistics from simulations
    # df_stat_sims_all_list <- vector("list", length = length(rows) * length(cols))
    # ind <- 0
    # start_time <- Sys.time()
    # for (row in rows) {
    #     for (col in cols) {
    #         ind <- ind + 1
    #         cat(paste("\nSTATISTICS FOR BATCH ", ind, "/", length(rows) * length(cols), "...\n", sep = ""))
    #         if (var1_name == "delta_selection") {
    #             filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_", var2_name, "=", scientific(var2_vals[col]))
    #             var1 <- scientific(var1_vals[1, row])
    #             var2 <- scientific(var2_vals[col])
    #         } else if (var2_name == "delta_selection") {
    #             filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
    #             var1 <- scientific(var1_vals[row])
    #             var2 <- scientific(var2_vals[1, col])
    #         } else {
    #             filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[col]))
    #             var1 <- scientific(var1_vals[row])
    #             var2 <- scientific(var2_vals[col])
    #         }
    #         if (compute_parallel == FALSE) {
    #             #--Get statistics for each simulation in sequential mode
    #             library(scales)
    #             library(vegan)
    #             df_stat_sims_list <- vector("list", n_simulations)
    #             for (sim in 1:n_simulations) {
    #                 filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
    #                 df_stat_sim <- statistics_multivar_one_simulation(filename, var1, var2, sim, plot_WGD)
    #                 gc()
    #                 df_stat_sims_list[[sim]] <- df_stat_sim
    #             }
    #         } else {
    #             #----Get statistics for each simulation in parallel mode
    #             library(parallel)
    #             library(pbapply)
    #             start_time <- Sys.time()
    #             #   Start parallel cluster
    #             if (is.null(n_cores)) {
    #                 numCores <- detectCores()
    #             } else {
    #                 numCores <- n_cores
    #             }
    #             cl <- makePSOCKcluster(numCores - 1)
    #             if (is.null(R_libPaths) == FALSE) {
    #                 R_libPaths <<- R_libPaths
    #                 clusterExport(cl, varlist = c(
    #                     "R_libPaths"
    #                 ))
    #                 clusterEvalQ(cl = cl, .libPaths(R_libPaths))
    #             }
    #             clusterEvalQ(cl = cl, library(scales))
    #             clusterEvalQ(cl = cl, library(vegan))
    #             #   Prepare input parameters for plotting
    #             var1 <<- var1
    #             var2 <<- var2
    #             statistics_multivar_one_simulation <<- statistics_multivar_one_simulation
    #             plot_WGD <<- plot_WGD
    #             compute_ploidy <<- compute_ploidy
    #             clusterExport(cl, varlist = c(
    #                 "filename_prefix",
    #                 "var1",
    #                 "var2",
    #                 "plot_WGD",
    #                 "statistics_multivar_one_simulation",
    #                 "compute_ploidy"
    #             ))
    #             #   Get statistics in parallel
    #             pbo <- pboptions(type = "txt")
    #             df_stat_sims_list <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(sim) {
    #                 filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
    #                 df_stat_sim <- statistics_multivar_one_simulation(filename, var1, var2, sim, plot_WGD)
    #                 return(df_stat_sim)
    #             })
    #             #   Stop parallel cluster
    #             stopCluster(cl)
    #             end_time <- Sys.time()
    #             print(end_time - start_time)
    #         }
    #         df_stat_sims_all_list[[ind]] <- rbindlist(df_stat_sims_list)
    #         # print(df_stat_sims_all_list[[ind]])
    #     }
    # }
    # end_time <- Sys.time()
    # print(end_time - start_time)
    # df_stat_sims_all <- rbindlist(df_stat_sims_all_list)
    # save(df_stat_sims_all, file = paste0(folder_workplace, "/", model_prefix, "_", "simulation_stats.rda"))
    load(paste0(folder_workplace, "/", model_prefix, "_", "simulation_stats.rda"))







    #-----------------------------------------Compute average statistics
    df_ploidy_dist <- data.frame(matrix(ncol = 7, nrow = 0))
    colnames(df_ploidy_dist) <- c("var1", "var2", "haploid", "diploid", "triploid", "tetraploid", "other")
    if (plot_WGD == TRUE) {
        df_WGD_clonality_dist <- data.frame(matrix(ncol = 5, nrow = 0))
        colnames(df_WGD_clonality_dist) <- c("var1", "var2", "clonal_WGD", "subclonal_WGD", "other")
    }

    df_stat_average <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df_stat_average) <- c("var1", "var2", "stat", "val")
    for (row in rows) {
        for (col in cols) {
            if (var1_name == "delta_selection") {
                filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_", var2_name, "=", scientific(var2_vals[col]))
                var1 <- scientific(var1_vals[1, row])
                var2 <- scientific(var2_vals[col])
            } else if (var2_name == "delta_selection") {
                filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
                var1 <- scientific(var1_vals[row])
                var2 <- scientific(var2_vals[1, col])
            } else {
                filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_", var1_name, "=", scientific(var1_vals[row]), "_", var2_name, "=", scientific(var2_vals[col]))
                var1 <- scientific(var1_vals[row])
                var2 <- scientific(var2_vals[col])
            }
            #---Statistics: clone count
            mean_clone_count <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clone_count")]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "clone_count", mean_clone_count)
            #---Statistics: fitness with respect to diploid clone
            mean_fitness <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "fitness")]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "fitness", mean_fitness)
            #---Statistics: Shannon diversity index
            mean_Shannon_index <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "Shannon_index")]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "Shannon_index", mean_Shannon_index)
            #---Statistics: ploidy
            mean_ploidy <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy")]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "ploidy", mean_ploidy)
            #---Statistics: distribution of ploidy
            simulations_ploidy <- round(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy")]))
            perc_1 <- 100 * length(which(simulations_ploidy == 1)) / length(simulations_ploidy)
            perc_2 <- 100 * length(which(simulations_ploidy == 2)) / length(simulations_ploidy)
            perc_3 <- 100 * length(which(simulations_ploidy == 3)) / length(simulations_ploidy)
            perc_4 <- 100 * length(which(simulations_ploidy == 4)) / length(simulations_ploidy)
            perc_other <- 100 * length(which(!(simulations_ploidy %in% c(1, 2, 3, 4)) == TRUE)) / length(simulations_ploidy)
            df_ploidy_dist[nrow(df_ploidy_dist) + 1, ] <- c(var1, var2, perc_1, perc_2, perc_3, perc_4, perc_other)

            #---Statistics: distribution of WGD clonality & subclonality
            if (plot_WGD == TRUE) {
                flag_clonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clonal_WGD")])
                flag_subclonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "subclonal_WGD")])
                flag_other <- as.numeric(flag_clonal_WGD == 0 & flag_subclonal_WGD == 0)
                perc_clonal_WGD <- 100 * sum(flag_clonal_WGD) / length(flag_clonal_WGD)
                perc_subclonal_WGD <- 100 * sum(flag_subclonal_WGD) / length(flag_subclonal_WGD)
                perc_other <- 100 * sum(flag_other) / length(flag_other)
                df_WGD_clonality_dist[nrow(df_WGD_clonality_dist) + 1, ] <- c(var1, var2, perc_clonal_WGD, perc_subclonal_WGD, perc_other)
            }
            # #---Statistics: percentage of clonal WGD
            # if (plot_WGD == TRUE) {
            #     vec_status_clonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clonal_WGD")])
            #     mean_percentage_clonal_WGD <- 100 * sum(vec_status_clonal_WGD) / length(vec_status_clonal_WGD)
            #     df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "clonal_WGD", mean_percentage_clonal_WGD)
            # }
            # #---Statistics: percentage of subclonal WGD
            # if (plot_WGD == TRUE) {
            #     vec_status_subclonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "subclonal_WGD")])
            #     mean_percentage_subclonal_WGD <- 100 * sum(vec_status_subclonal_WGD) / length(vec_status_subclonal_WGD)
            #     df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "subclonal_WGD", mean_percentage_subclonal_WGD)
            # }
            # #---Statistics: percentage of WGD (both clonal and subclonal)
            # if (plot_WGD == TRUE) {
            #     vec_status_clonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clonal_WGD")])
            #     vec_status_subclonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "subclonal_WGD")])
            #     vec_status_WGD <- vec_status_clonal_WGD + vec_status_subclonal_WGD
            #     vec_status_WGD[which(vec_status_WGD > 1)] <- 1
            #     mean_percentage_WGD <- 100 * sum(vec_status_WGD) / length(vec_status_WGD)
            #     df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "WGD", mean_percentage_WGD)
            # }
            #---Statistics: event count before clonal WGD
            if (plot_WGD == TRUE) {
                vec_event_count_pre_clonal_WGD <- as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clonal_WGD")])
                percent_NA <- length(is.na(vec_event_count_pre_clonal_WGD)) / length(vec_event_count_pre_clonal_WGD)
                if (percent_NA < 0.1) {
                    mean_event_count_pre_clonal_WGD <- NA
                } else {
                    vec_event_count_pre_clonal_WGD <- vec_event_count_pre_clonal_WGD[which(!is.na(vec_event_count_pre_clonal_WGD))]
                    mean_event_count_pre_clonal_WGD <- 100 * sum(vec_event_count_pre_clonal_WGD) / length(vec_event_count_pre_clonal_WGD)
                }
                df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "event_count_pre_clonal_WGD", mean_event_count_pre_clonal_WGD)
            }
        }
    }
    save(df_stat_average, file = paste0(folder_workplace, "/", model_prefix, "_", "average_stats.rda"))
    load(paste0(folder_workplace, "/", model_prefix, "_", "average_stats.rda"))
    save(df_ploidy_dist, file = paste0(folder_workplace, "/", model_prefix, "_", "ploidy_distribution.rda"))
    load(paste0(folder_workplace, "/", model_prefix, "_", "ploidy_distribution.rda"))
    if (plot_WGD == TRUE) {
        save(df_WGD_clonality_dist, file = paste0(folder_workplace, "/", model_prefix, "_", "WGD_clonality_distribution.rda"))
        load(paste0(folder_workplace, "/", model_prefix, "_", "ploidy_distribution.rda"))
    }







    #-------------------------------------Prepare dataframe for plotting
    df_stat_average$val <- as.numeric(df_stat_average$val)

    df_ploidy_dist$var1 <- as.numeric(df_ploidy_dist$var1)
    df_ploidy_dist$var2 <- as.numeric(df_ploidy_dist$var2)
    df_ploidy_dist$haploid <- as.numeric(df_ploidy_dist$haploid)
    df_ploidy_dist$diploid <- as.numeric(df_ploidy_dist$diploid)
    df_ploidy_dist$triploid <- as.numeric(df_ploidy_dist$triploid)
    df_ploidy_dist$tetraploid <- as.numeric(df_ploidy_dist$tetraploid)
    df_ploidy_dist$other <- as.numeric(df_ploidy_dist$other)

    if (plot_WGD == TRUE) {
        df_WGD_clonality_dist$var1 <- as.numeric(df_WGD_clonality_dist$var1)
        df_WGD_clonality_dist$var2 <- as.numeric(df_WGD_clonality_dist$var2)
        df_WGD_clonality_dist$clonal_WGD <- as.numeric(df_WGD_clonality_dist$clonal_WGD)
        df_WGD_clonality_dist$subclonal_WGD <- as.numeric(df_WGD_clonality_dist$subclonal_WGD)
        df_WGD_clonality_dist$other <- as.numeric(df_WGD_clonality_dist$other)
    }
    #   Plot fitness on logscale
    df_stat_average$val[which(df_stat_average$stat == "fitness")] <- log(df_stat_average$val[which(df_stat_average$stat == "fitness")])
    #   Sort dataframe according to increasing variables
    if (is.vector(var1_vals)) {
        df_stat_average$var1 <- factor(scientific(as.numeric(df_stat_average$var1)), levels = scientific(sort(var1_vals)))
        x_ticks_breaks <- sort(var1_vals)
    } else if (is.matrix(var1_vals)) {
        df_stat_average$var1 <- factor(scientific(as.numeric(df_stat_average$var1)), levels = scientific(sort(var1_vals[1, ])))
        x_ticks_breaks <- sort(var1_vals[1, ])
    }
    if (is.vector(var2_vals)) {
        df_stat_average$var2 <- factor(scientific(as.numeric(df_stat_average$var2)), levels = scientific(sort(var2_vals)))
        y_ticks_breaks <- sort(var2_vals)
    } else if (is.matrix(var2_vals)) {
        df_stat_average$var2 <- factor(scientific(as.numeric(df_stat_average$var2)), levels = scientific(sort(var2_vals[1, ])))
        y_ticks_breaks <- sort(var2_vals[1, ])
    }
    df_stat_average <- df_stat_average[which((is.na(df_stat_average$var1) == FALSE) & (is.na(df_stat_average$var2) == FALSE)), ]
    #   Force identical scale for distribution maps
    x_ticks_breaks <- sort(unique(df_ploidy_dist$var1))
    y_ticks_breaks <- sort(unique(df_ploidy_dist$var2))
    for (ind in 1:length(x_ticks_breaks)) {
        df_ploidy_dist$var1[which(df_ploidy_dist$var1 == x_ticks_breaks[ind])] <- ind
        if (plot_WGD == TRUE) df_WGD_clonality_dist$var1[which(df_WGD_clonality_dist$var1 == x_ticks_breaks[ind])] <- ind
    }
    for (ind in 1:length(y_ticks_breaks)) {
        df_ploidy_dist$var2[which(df_ploidy_dist$var2 == y_ticks_breaks[ind])] <- ind
        if (plot_WGD == TRUE) df_WGD_clonality_dist$var2[which(df_WGD_clonality_dist$var2 == y_ticks_breaks[ind])] <- ind
    }

    #---------------------------------------Plot statistics: clone count
    filename <- paste0(model_prefix, "_1_clone_count.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "clone_count"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_viridis(discrete = FALSE, name = "Clone count") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #---------------------------Plot statistics: Shannon diversity index
    filename <- paste0(model_prefix, "_2_Shannon_index.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "Shannon_index"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_viridis(discrete = FALSE, name = "Shannon index") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #-------------Plot statistics: fitness with respect to diploid clone
    filename <- paste0(model_prefix, "_3_fitness.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "fitness"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "RdPu", name = "log(Fitness)") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #--------------------------------------------Plot statistics: ploidy
    filename <- paste0(model_prefix, "_4_ploidy.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "ploidy"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_gradientn(
            colours = c("black", "blue3", "white", "forestgreen", "firebrick2", "red4"),
            breaks = c(0, 1, 2, 3, 4, 5),
            labels = c(0, 1, 2, 3, 4, 5),
            limits = c(0, 5),
            name = "ploidy"
        ) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()



    #----------------------------Plot statistics: distribution of ploidy
    filename <- paste0(model_prefix, "_4_ploidy_distribution.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    p <- ggplot() +
        geom_scatterpie(aes(x = var1, y = var2, r = 0.45),
            data = df_ploidy_dist,
            cols = c("haploid", "diploid", "triploid", "tetraploid", "other")
        ) +
        coord_equal() +
        labs(fill = "percentage") +
        xlab(var1_lab) +
        scale_x_continuous(breaks = 1:length(x_ticks_breaks), labels = scientific(x_ticks_breaks)) +
        ylab(var2_lab) +
        scale_y_continuous(breaks = 1:length(y_ticks_breaks), labels = scientific(y_ticks_breaks)) +
        scale_fill_manual(values = c("blue3", "white", "forestgreen", "firebrick2", "black")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #----------------Plot statistics: distribution of WGD (sub)clonality
    if (plot_WGD == TRUE) {
        filename <- paste0(model_prefix, "_5_WGD_clonality_distribution.jpeg")
        jpeg(file = filename, width = 1000, height = 1100)
        # df_WGD_clonality_dist$radius<-
        p <- ggplot() +
            geom_scatterpie(aes(x = var1, y = var2, r = 0.45),
                data = df_WGD_clonality_dist,
                cols = c("clonal_WGD", "subclonal_WGD")
            ) +
            coord_equal() +
            labs(fill = "percentage") +
            xlab(var1_lab) +
            scale_x_continuous(breaks = 1:length(x_ticks_breaks), labels = scientific(x_ticks_breaks)) +
            ylab(var2_lab) +
            scale_y_continuous(breaks = 1:length(y_ticks_breaks), labels = scientific(y_ticks_breaks)) +
            scale_fill_manual(values = c("sienna", "wheat")) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
        print(p)
        dev.off()
    }



    # #--------------------------Plot statistics: percentage of clonal WGD
    # if (plot_WGD == TRUE) {
    #     filename <- paste0(model_prefix, "_5_clonal_WGD.jpeg")
    #     jpeg(file = filename, width = 1000, height = 1100)
    #     df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "clonal_WGD"), ]
    #     p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
    #         geom_tile() +
    #         coord_equal() +
    #         xlab(var1_lab) +
    #         ylab(var2_lab) +
    #         scale_fill_gradientn(
    #             colours = c("black", "sienna", "goldenrod1", "wheat", "white"),
    #             breaks = c(0, 25, 50, 75, 100),
    #             labels = c(0, 25, 50, 75, 100),
    #             limits = c(0, 100),
    #             name = "Percentage of clonal WGD"
    #         ) +
    #         theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    #     print(p)
    #     dev.off()
    # }
    # #-----------------------Plot statistics: percentage of subclonal WGD
    # if (plot_WGD == TRUE) {
    #     filename <- paste0(model_prefix, "_6_subclonal_WGD.jpeg")
    #     jpeg(file = filename, width = 1000, height = 1100)
    #     df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "subclonal_WGD"), ]
    #     p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
    #         geom_tile() +
    #         coord_equal() +
    #         xlab(var1_lab) +
    #         ylab(var2_lab) +
    #         scale_fill_gradientn(
    #             colours = c("black", "sienna", "goldenrod1", "wheat", "white"),
    #             breaks = c(0, 25, 50, 75, 100),
    #             labels = c(0, 25, 50, 75, 100),
    #             limits = c(0, 100),
    #             name = "Percentage of subclonal WGD"
    #         ) +
    #         theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    #     print(p)
    #     dev.off()
    # }
    # #---------------------------------Plot statistics: percentage of WGD
    # if (plot_WGD == TRUE) {
    #     filename <- paste0(model_prefix, "_6_WGD.jpeg")
    #     jpeg(file = filename, width = 1000, height = 1100)
    #     df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "WGD"), ]
    #     p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
    #         geom_tile() +
    #         coord_equal() +
    #         xlab(var1_lab) +
    #         ylab(var2_lab) +
    #         scale_fill_gradientn(
    #             colours = c("black", "sienna", "goldenrod1", "wheat", "white"),
    #             breaks = c(0, 25, 50, 75, 100),
    #             labels = c(0, 25, 50, 75, 100),
    #             limits = c(0, 100),
    #             name = "Percentage of WGD (clonal or subclonal)"
    #         ) +
    #         theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    #     print(p)
    #     dev.off()
    # }
    #---------------------Plot statistics: event count before clonal WGD
    if (plot_WGD == TRUE) {
        filename <- paste0(model_prefix, "_6_event_count_pre_clonal_WGD.jpeg")
        jpeg(file = filename, width = 1000, height = 1100)
        df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "event_count_pre_clonal_WGD"), ]
        p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
            geom_tile() +
            coord_equal() +
            xlab(var1_lab) +
            ylab(var2_lab) +
            scale_fill_gradientn(
                colours = c("white", "chartreuse4"),
                na.value = "black",
                name = "Event count before clonal WGD"
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
        print(p)
        dev.off()
    }
}

statistics_multivar_one_simulation <- function(filename, var1, var2, sim, plot_WGD) {
    load(filename)
    #--------------------------------Create dataframe for all statistics
    df_stat_sim <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_stat_sim) <- c("var1", "var2", "sim", "stat", "val")
    #---------------------------Input required variables from simulation
    genotype_list_selection_rate <- simulation$clonal_evolution$genotype_list_selection_rate
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes



    table_clone <- as.data.frame(table(ID = simulation$sample$sample_clone_ID))
    Clone_ID <- as.numeric(as.vector(table_clone$ID))
    #--------------------------------------------Statistics: clone count
    clone_count <- length(Clone_ID)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clone_count", clone_count)
    #------------------Statistics: fitness with respect to diploid clone
    table_clone$fitness <- genotype_list_selection_rate[Clone_ID] / genotype_list_selection_rate[1]
    fitness <- sum((table_clone$fitness) * (table_clone$Freq)) / sum(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "fitness", fitness)
    #--------------------------------Statistics: Shannon diversity index
    Shannon_index <- diversity(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "Shannon_index", Shannon_index)
    #-------------Statistics: CCF for diploid+5 biggest subclones+others







    #-------------------------------------------------Statistics: ploidy
    N_chromosomes <- length(genotype_list_ploidy_chrom[[1]])
    vec_CN_block_no <- rep(0, N_chromosomes)
    for (chrom in 1:N_chromosomes) {
        vec_CN_block_no[chrom] <- length(genotype_list_ploidy_block[[1]][[chrom]][[1]])
    }
    #   Compute ploidy for each unique clone
    table_clone$ploidy <- 0
    for (i in 1:nrow(table_clone)) {
        table_clone$ploidy[i] <- compute_ploidy(vec_CN_block_no, genotype_list_ploidy_chrom[[Clone_ID[i]]], genotype_list_ploidy_block[[Clone_ID[i]]])
    }
    #   Compute mean ploidy for simulation
    ploidy <- sum((table_clone$ploidy) * (table_clone$Freq)) / sum(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "ploidy", ploidy)
    #-----------------Statistics: classification as clonal/subclonal WGD
    if (plot_WGD == TRUE) {
        table_clone$WGD <- 0
        #   Find WGD status of each unique clone
        for (i in 1:nrow(table_clone)) {
            clone <- Clone_ID[i]
            WGD_status <- 0
            clone_node <- clone
            while (clone_node > 0) {
                events <- evolution_genotype_changes[[clone_node]]
                if (length(events) > 0) {
                    for (event in 1:length(events)) {
                        if (events[[event]][1] == "whole-genome-duplication") WGD_status <- 1
                    }
                }
                if (WGD_status == 1) break
                clone_node <- evolution_origin[clone_node]
            }
            table_clone$WGD[i] <- WGD_status
        }
        #   Define simulation as clonal WGD or subclonal WGD or neither
        if ((min(table_clone$WGD) == 0) & (max(table_clone$WGD) == 1)) {
            clonal_WGD <- 0
            subclonal_WGD <- 1
        } else if ((min(table_clone$WGD) == 1) & (max(table_clone$WGD) == 1)) {
            clonal_WGD <- 1
            subclonal_WGD <- 0
        } else {
            clonal_WGD <- 0
            subclonal_WGD <- 0
        }
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clonal_WGD", clonal_WGD)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "subclonal_WGD", subclonal_WGD)
    }
    #-----Statistics: event count before clonal WGD
    if (plot_WGD == TRUE) {
        if (clonal_WGD == 1) {
            table_clone$event_count_pre_first_WGD <- 0
            #   Find number of events before 1st WGD in each unique clone
            for (i in 1:nrow(table_clone)) {
                clone <- Clone_ID[i]
                #   Find clonal lineage
                clone_lineage <- c()
                clone_node <- clone
                while (clone_node > 0) {
                    clone_lineage <- c(clone_node, clone_lineage)
                    clone_node <- evolution_origin[clone_node]
                }
                clone_lineage <- sort(clone_lineage)
                #   Find number of events before 1st WGD
                event_count <- 0
                for (j in 1:length(clone_lineage)) {
                    clone_node <- clone_lineage[j]
                    events <- evolution_genotype_changes[[clone_node]]
                    WGD <- FALSE
                    if (length(events) > 0) {
                        for (event in 1:length(events)) {
                            if (events[[event]][1] == "whole-genome-duplication") {
                                WGD <- TRUE
                            } else if (WGD == FALSE) {
                                event_count <- event_count + 1
                            }
                        }
                    }
                    if (WGD == TRUE) break
                }
                table_clone$event_count_pre_first_WGD[i] <- event_count
            }
            event_count_before_WGD <- sum((table_clone$event_count_pre_first_WGD) * (table_clone$Freq)) / sum(table_clone$Freq)
        } else {
            event_count_before_WGD <- NA
        }
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "event_count_pre_clonal_WGD", event_count_before_WGD)
    }





    #--------------------------Return the statistics for this simulation
    return(df_stat_sim)
}

compute_ploidy <- function(vec_CN_block_no, ploidy_chrom, ploidy_block) {
    N_chromosomes <- length(vec_CN_block_no)
    vec_CN_all <- c()
    for (chrom in 1:N_chromosomes) {
        vec_CN <- rep(0, vec_CN_block_no[chrom])
        no_strands <- ploidy_chrom[chrom]
        if (no_strands > 0) {
            for (strand in 1:no_strands) {
                vec_CN <- vec_CN + ploidy_block[[chrom]][[strand]]
            }
        }
        vec_CN_all <- c(vec_CN_all, vec_CN)
    }
    ploidy <- mean(vec_CN_all)
    return(ploidy)
}

#' @export
simulation_statistics_multivar <- function(model_prefix = "",
                                           folder_workplace = NULL,
                                           folder_plots = NULL,
                                           model_variables_base,
                                           var1_name = "",
                                           var1_vals = c(),
                                           var2_name = "",
                                           var2_vals = c(),
                                           n_simulations = 0,
                                           compute_parallel = TRUE,
                                           n_cores = NULL) {
    if (is.vector(var1_vals)) {
        inds <- 1:length(var1_vals)
    } else if (is.matrix(var1_vals)) {
        inds <- 1:ncol(var1_vals)
    }
    #------------------------------------Get statistics from simulations
    plot_name <- c("A", "B", "C", "D")

    for (ind in inds) {
        cat(paste("\nSTATISTICS FOR BATCH ", ind, "/", length(inds), "...\n", sep = ""))
        if ((var1_name == "delta_selection") & (is.matrix(var1_vals))) {
            filename_prefix <<- paste0(model_prefix, "_", var1_name, "=", scientific(var1_vals[1, ind]), "&", scientific(var1_vals[2, ind]), "_", var2_name, "=", scientific(var2_vals[ind]))
        } else if ((var1_name == "delta_selection") & (is.vector(var1_vals))) {
            filename_prefix <<- paste0(model_prefix, "_", var1_name, "=", scientific(var1_vals[1]), "&", scientific(var1_vals[2]), "_", var2_name, "=", scientific(var2_vals[ind]))
        } else if ((var2_name == "delta_selection") & (is.matrix(var2_vals))) {
            filename_prefix <<- paste0(model_prefix, "_", var1_name, "=", scientific(var1_vals[ind]), "_", var2_name, "=", scientific(var2_vals[1, ind]), "&", scientific(var2_vals[2, ind]))
        } else if ((var2_name == "delta_selection") & (is.vector(var2_vals))) {
            filename_prefix <<- paste0(model_prefix, "_", var1_name, "=", scientific(var1_vals[ind]), "_", var2_name, "=", scientific(var2_vals[1]), "&", scientific(var2_vals[2]))
        } else {
            filename_prefix <<- paste0(model_prefix, "_", var1_name, "=", scientific(var1_vals[ind]), "_", var2_name, "=", scientific(var2_vals[ind]))
        }
        print(filename_prefix)
        #--------------------------Make plots for individual simulations
        plot_all(
            model = filename_prefix,
            n_simulations = n_simulations,
            unit_time = "day",
            folder_workplace = folder_workplace,
            folder_plots = folder_plots,
            compute_parallel = compute_parallel
        )
    }
}
