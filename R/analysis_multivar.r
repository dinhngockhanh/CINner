#' @export
simulator_multivar <- function(model_variables_base = list(),
                               model_prefix = "",
                               var1_name = "",
                               var1_vals = c(),
                               var1_labs = NULL,
                               var2_name = "",
                               var2_vals = c(),
                               var2_labs = NULL,
                               extra_var = NULL,
                               n_simulations = 0,
                               stage_final = 3,
                               n_clones_min = 0,
                               n_clones_max = Inf,
                               save_simulation = TRUE,
                               lite_memory = TRUE,
                               neutral_variations = FALSE,
                               internal_nodes_cn_info = FALSE,
                               save_newick_tree = FALSE,
                               save_cn_profile = FALSE,
                               save_cn_clones = FALSE,
                               build_cn = FALSE,
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
            if (is.vector(var1_vals)) {
                var1_val <- var1_vals[row]
            } else if (is.matrix(var1_vals)) {
                var1_val <- var1_vals[, row]
            }
            if (is.vector(var2_vals)) {
                var2_val <- var2_vals[col]
            } else if (is.matrix(var2_vals)) {
                var2_val <- var2_vals[, col]
            }
            #-----------------------Create model variables for the batch
            #   Initialize parameter set
            model_variables <- model_variables_base
            #   Fix variable 1 in parameter set
            if (var1_name %in% model_variables$general_variables$Variable) {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == var1_name)] <- var1_val
            } else if (var1_name %in% model_variables$selection_model$Variable) {
                model_variables$selection_model$Value[which(model_variables$selection_model$Variable == var1_name)] <- var1_val
            } else if (var1_name == "scale_selection") {
                s_rate_scale <- var1_val
                s_rate_base <- model_variables$chromosome_arm_library$s_rate
                s_rate <- s_rate_base
                s_rate[which(s_rate_base > 1)] <- 1 + s_rate_scale * (s_rate_base[which(s_rate_base > 1)] - 1)
                s_rate[which(s_rate_base < 1)] <- 1 / (1 + s_rate_scale * (1 / s_rate_base[which(s_rate_base < 1)] - 1))
                model_variables$chromosome_arm_library$s_rate <- s_rate
            } else if (var1_name == "scale_selection_gain") {
                s_rate_scale <- var1_val
                s_rate_base <- model_variables$chromosome_arm_library$s_rate
                locs <- which(s_rate_base > 1)
                model_variables$chromosome_arm_library$s_rate[locs] <- 1 + s_rate_scale * (s_rate_base[locs] - 1)
            } else if (var1_name == "scale_selection_loss") {
                s_rate_scale <- var1_val
                s_rate_base <- model_variables$chromosome_arm_library$s_rate
                locs <- which(s_rate_base < 1)
                model_variables$chromosome_arm_library$s_rate[locs] <- 1 / (1 + s_rate_scale * (1 / s_rate_base[locs] - 1))
            } else if (var1_name == "delta_selection") {
                delta_sel_gain_per_arm <- sqrt(1 + var1_val[1]) - 1
                delta_sel_loss_per_arm <- sqrt(1 + var1_val[2]) - 1
                chrom_gains <- sample(model_variables_base$cn_info$Chromosome, round(length(model_variables_base$cn_info$Chromosome) / 2), replace = FALSE)
                chrom_losses <- setdiff(model_variables_base$cn_info$Chromosome, chrom_gains)
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_gains)] <- 1 + delta_sel_gain_per_arm
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_losses)] <- 1 / (1 + delta_sel_loss_per_arm)
            } else if (var1_name == "vec_cell_count") {
                model_variables$population_dynamics$Total_cell_count <- var1_val * sum(model_variables$population_dynamics$Total_cell_count) / sum(var1_val)
            } else if (var1_name == "scale_cell_count") {
                model_variables$population_dynamics$Total_cell_count <- model_variables$population_dynamics$Total_cell_count * var1_val * length(model_variables$population_dynamics$Total_cell_count) / sum(model_variables$population_dynamics$Total_cell_count)
            } else {
                stop("VARIABLE 1 IS NOT RECOGNIZED")
            }
            #   Fix variable 2 in parameter set
            if (var2_name %in% model_variables$general_variables$Variable) {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == var2_name)] <- var2_val
            } else if (var2_name %in% model_variables$selection_model$Variable) {
                model_variables$selection_model$Value[which(model_variables$selection_model$Variable == var2_name)] <- var2_val
            } else if (var2_name == "scale_selection") {
                s_rate_base <- model_variables_base$chromosome_arm_library$s_rate
                s_rate_scale <- var2_val
                s_rate <- s_rate_base
                s_rate[which(s_rate_base > 1)] <- 1 + s_rate_scale * (s_rate_base[which(s_rate_base > 1)] - 1)
                s_rate[which(s_rate_base < 1)] <- 1 / (1 + s_rate_scale * (1 / s_rate_base[which(s_rate_base < 1)] - 1))
                model_variables$chromosome_arm_library$s_rate <- s_rate
            } else if (var2_name == "scale_selection_gain") {
                s_rate_scale <- var2_val
                s_rate_base <- model_variables$chromosome_arm_library$s_rate
                locs <- which(s_rate_base > 1)
                model_variables$chromosome_arm_library$s_rate[locs] <- 1 + s_rate_scale * (s_rate_base[locs] - 1)
            } else if (var2_name == "scale_selection_loss") {
                s_rate_scale <- var2_val
                s_rate_base <- model_variables$chromosome_arm_library$s_rate
                locs <- which(s_rate_base < 1)
                model_variables$chromosome_arm_library$s_rate[locs] <- 1 / (1 + s_rate_scale * (1 / s_rate_base[locs] - 1))
            } else if (var2_name == "delta_selection") {
                delta_sel_gain_per_arm <- sqrt(1 + var2_val[1]) - 1
                delta_sel_loss_per_arm <- sqrt(1 + var2_val[2]) - 1
                chrom_gains <- sample(model_variables_base$cn_info$Chromosome, round(length(model_variables_base$cn_info$Chromosome) / 2), replace = FALSE)
                chrom_losses <- setdiff(model_variables_base$cn_info$Chromosome, chrom_gains)
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_gains)] <- 1 + delta_sel_gain_per_arm
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Chromosome %in% chrom_losses)] <- 1 / (1 + delta_sel_loss_per_arm)
            } else if (var2_name == "vec_cell_count") {
                model_variables$population_dynamics$Total_cell_count <- var2_val * sum(model_variables$population_dynamics$Total_cell_count) / sum(var2_val)
            } else if (var2_name == "scale_cell_count") {
                model_variables$population_dynamics$Total_cell_count <- model_variables$population_dynamics$Total_cell_count * var2_val * length(model_variables$population_dynamics$Total_cell_count) / sum(model_variables$population_dynamics$Total_cell_count)
            } else {
                stop("VARIABLE 2 IS NOT RECOGNIZED")
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
            model_variables <- CHECK_model_variables(model_variables)
            #   Save model variables
            model_name <- paste0(model_prefix, "_")
            if (!is.null(var1_labs)) {
                model_name <- paste0(model_name, var1_name, "=", var1_labs[row], "_")
            } else if (is.vector(var1_vals)) {
                model_name <- paste0(model_name, var1_name, "=", scientific(var1_vals[row]), "_")
            } else if (is.matrix(var1_vals)) {
                model_name <- paste0(model_name, var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_")
            }
            if (!is.null(var2_labs)) {
                model_name <- paste0(model_name, var2_name, "=", var2_labs[col])
            } else if (is.vector(var2_vals)) {
                model_name <- paste0(model_name, var2_name, "=", scientific(var2_vals[col]))
            } else if (is.matrix(var2_vals)) {
                model_name <- paste0(model_name, var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
            }
            SAVE_model_variables(model_name = model_name, model_variables = model_variables)
            #   Create simulations
            cat("=======================================================\n")
            cat("=======================================================\n")
            cat("=======================================================\n")
            cat(paste("\nSIMULATIONS FOR BATCH ", ind, "/", length(rows) * length(cols), "...\n", sep = ""))
            if (is.matrix(var1_vals)) {
                cat(paste(var1_name, " = ", scientific(var1_vals[1, row]), " & ", scientific(var1_vals[2, row]), "\n", sep = ""))
            } else {
                cat(paste(var1_name, " = ", scientific(var1_vals[row]), "\n", sep = ""))
            }
            if (is.matrix(var2_vals)) {
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
                build_cn = build_cn,
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
                lite_memory = lite_memory,
                R_libPaths = R_libPaths
            )
            end_time <- Sys.time()
            print(end_time - start_time)
            cat("\n")
            #-----------------------------Plot simulations for the batch
            if (plot == TRUE) {
                start_time <- Sys.time()
                # plot_all(
                #     model = model_name,
                #     n_simulations = n_simulations,
                #     unit_time = model_variables$general_variables$Unit[which(model_variables$general_variables$Variable == "T_end_time")],
                #     folder_workplace = folder_workplace,
                #     folder_plots = folder_workplace,
                #     compute_parallel = compute_parallel,
                #     R_libPaths = R_libPaths
                # )

                plot_cn_heatmap(
                    model = model_name,
                    n_simulations = n_simulations,
                    folder_workplace = folder_workplace,
                    folder_plots = folder_workplace,
                    plotcol = "total-copy",
                    CN_data = "TRUTH",
                    phylo = "TRUTH",
                    width = 1000,
                    height = 1000,
                    compute_parallel = compute_parallel,
                    n_cores = n_cores,
                    R_libPaths = R_libPaths
                )

                # plot_cn_heatmap(
                #     model = model_name,
                #     n_simulations = n_simulations,
                #     folder_workplace = folder_workplace,
                #     folder_plots = folder_workplace,
                #     plotcol = "total-copy",
                #     CN_data = "NEUTRAL-VARIATIONS",
                #     phylo = "TRUTH",
                #     width = 1000,
                #     height = 1000,
                #     compute_parallel = compute_parallel,
                #     n_cores = n_cores,
                #     R_libPaths = R_libPaths
                # )

                end_time <- Sys.time()
                print(end_time - start_time)
                cat("\n")
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
statistics_multivar_matrix <- function(model_prefix = "",
                                       folder_workplace = "",
                                       var1_name = "",
                                       var1_vals = c(),
                                       var1_labs = NULL,
                                       var2_name = "",
                                       var2_vals = c(),
                                       var2_labs = NULL,
                                       n_simulations = 0,
                                       plot_WGD = FALSE,
                                       plot_misseg = FALSE,
                                       compute_parallel = TRUE,
                                       n_cores = NULL,
                                       R_libPaths = NULL) {
    library(data.table)
    library(scales)
    library(vegan)
    library(ggplot2)
    library(viridis)
    library(hrbrthemes)
    library(scatterpie)
    if (var1_name == "prob_CN_whole_genome_duplication") {
        var1_lab <- "Probability of WGD"
    } else if (var1_name == "alpha_aneuploidy") {
        var1_lab <- "WGD-aneuploidy rate"
    } else if (var1_name == "scale_selection") {
        var1_lab <- "Scale of selection rates"
    } else if (var1_name == "scale_selection_gain") {
        var1_lab <- "Scale of selection rates for ONCOGENE arms"
    } else if (var1_name == "scale_selection_loss") {
        var1_lab <- "Scale of selection rates for TSG arms"
    } else if (var1_name == "delta_selection") {
        var1_lab <- "Selection rate"
    } else if (var1_name == "prob_CN_missegregation") {
        var1_lab <- "Probability of missegregation"
    } else if (var1_name == "vec_cell_count") {
        var1_lab <- "Growth mode"
    } else if (var1_name == "scale_cell_count") {
        var1_lab <- "Average cell count"
    } else {
        stop("VARIABLE 1 IS NOT RECOGNIZED")
    }
    if (var2_name == "prob_CN_whole_genome_duplication") {
        var2_lab <- "Probability of WGD"
    } else if (var2_name == "alpha_aneuploidy") {
        var2_lab <- "WGD-aneuploidy rate"
    } else if (var2_name == "scale_selection") {
        var2_lab <- "Scale of selection rates"
    } else if (var2_name == "scale_selection_gain") {
        var2_lab <- "Scale of selection rates for ONCOGENE arms"
    } else if (var2_name == "scale_selection_loss") {
        var2_lab <- "Scale of selection rates for TSG arms"
    } else if (var2_name == "delta_selection") {
        var2_lab <- "Selection rate"
    } else if (var2_name == "prob_CN_missegregation") {
        var2_lab <- "Probability of missegregation"
    } else if (var2_name == "vec_cell_count") {
        var2_lab <- "Growth mode"
    } else if (var2_name == "scale_cell_count") {
        var2_lab <- "Average cell count"
    } else {
        stop("VARIABLE 2 IS NOT RECOGNIZED")
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
    #------------------------------------Get statistics from simulations
    ind <- 0
    start_time <- Sys.time()
    for (row in rows) {
        for (col in cols) {
            ind <- ind + 1
            cat(paste("\nSTATISTICS FOR BATCH ", ind, "/", length(rows) * length(cols), "...\n", sep = ""))
            filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_")
            if (!is.null(var1_labs)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", var1_labs[row], "_")
                var1 <- var1_labs[row]
            } else if (is.vector(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[row]), "_")
                var1 <- scientific(var1_vals[row])
            } else if (is.matrix(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_")
                var1 <- scientific(var1_vals[1, row])
            }
            if (!is.null(var2_labs)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", var2_labs[col])
                var2 <- var2_labs[col]
            } else if (is.vector(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[col]))
                var2 <- scientific(var2_vals[col])
            } else if (is.matrix(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
                var2 <- scientific(var2_vals[1, col])
            }
            if (compute_parallel == FALSE) {
                library(signals)
                #--Get statistics for each simulation in sequential mode
                df_stat_sims_list <- vector("list", n_simulations)
                for (sim in 1:n_simulations) {
                    filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
                    df_stat_sim <- statistics_multivar_one_simulation(filename, var1_name, var2_name, var1, var2, sim, plot_WGD, plot_misseg)
                    df_stat_sims_list[[sim]] <- df_stat_sim
                }
            } else {
                #----Get statistics for each simulation in parallel mode
                library(parallel)
                library(pbapply)
                start_time <- Sys.time()
                #   Start parallel cluster
                if (is.null(n_cores)) {
                    numCores <- detectCores()
                } else {
                    numCores <- n_cores
                }
                cl <- makePSOCKcluster(numCores - 1)
                if (is.null(R_libPaths) == FALSE) {
                    R_libPaths <<- R_libPaths
                    clusterExport(cl, varlist = c(
                        "R_libPaths"
                    ))
                    clusterEvalQ(cl = cl, .libPaths(R_libPaths))
                }
                clusterEvalQ(cl = cl, library(scales))
                clusterEvalQ(cl = cl, library(vegan))
                clusterEvalQ(cl = cl, library(signals))
                #   Prepare input parameters for plotting
                var1_name <<- var1_name
                var2_name <<- var2_name
                var1 <<- var1
                var2 <<- var2
                # statistics_multivar_one_simulation <<- statistics_multivar_one_simulation
                plot_WGD <<- plot_WGD
                plot_misseg <<- plot_misseg
                clusterExport(cl, varlist = c(
                    "filename_prefix",
                    "var1_name",
                    "var2_name",
                    "var1",
                    "var2",
                    "plot_WGD",
                    "statistics_multivar_one_simulation",
                    "get_cn_profile", "normalize_cell_ploidy", "calc_state_mode"
                ))
                #   Get statistics in parallel
                pbo <- pboptions(type = "txt")
                df_stat_sims_list <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(sim) {
                    filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
                    df_stat_sim <- statistics_multivar_one_simulation(filename, var1_name, var2_name, var1, var2, sim, plot_WGD, plot_misseg)
                    return(df_stat_sim)
                })
                #   Stop parallel cluster
                stopCluster(cl)
                end_time <- Sys.time()
                print(end_time - start_time)
            }
            df_stat_sims <- rbindlist(df_stat_sims_list)
            filename <- paste0(filename_prefix, "_simulation_stats.rda")
            save(df_stat_sims, file = filename)
        }
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    #---------Combine statistics from all simulations into one dataframe
    df_stat_sims_all_list <- vector("list", length = length(rows) * length(cols))
    ind <- 0
    for (row in rows) {
        for (col in cols) {
            ind <- ind + 1
            filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_")
            if (!is.null(var1_labs)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", var1_labs[row], "_")
                var1 <- var1_labs[row]
            } else if (is.vector(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[row]), "_")
                var1 <- scientific(var1_vals[row])
            } else if (is.matrix(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_")
                var1 <- scientific(var1_vals[1, row])
            }
            if (!is.null(var2_labs)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", var2_labs[col])
                var2 <- var2_labs[col]
            } else if (is.vector(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[col]))
                var2 <- scientific(var2_vals[col])
            } else if (is.matrix(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
                var2 <- scientific(var2_vals[1, col])
            }
            filename <- paste0(filename_prefix, "_simulation_stats.rda")
            load(filename)
            df_stat_sims_all_list[[ind]] <- df_stat_sims
        }
    }
    df_stat_sims_all <- rbindlist(df_stat_sims_all_list)
    save(df_stat_sims_all, file = paste0(folder_workplace, "/", model_prefix, "_", "simulation_stats.rda"))
    #-----------------------------------------Compute average statistics
    load(paste0(folder_workplace, "/", model_prefix, "_", "simulation_stats.rda"))
    df_ploidy_dist <- data.frame(matrix(ncol = 7, nrow = 0))
    colnames(df_ploidy_dist) <- c("var1", "var2", "haploid", "diploid", "triploid", "tetraploid", "other")
    df_nonviability <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_nonviability) <- c("var1", "var2", "diploid", "tetraploid", "other")
    if (plot_misseg == TRUE) {
        df_missegregation <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(df_missegregation) <- c("var1", "var2", "clonal_count_missegregation", "subclonal_count_missegregation")
    }
    if (plot_WGD == TRUE) {
        df_WGD_clonality_dist <- data.frame(matrix(ncol = 5, nrow = 0))
        colnames(df_WGD_clonality_dist) <- c("var1", "var2", "clonal_WGD", "subclonal_WGD", "other")
    }
    df_stat_average <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df_stat_average) <- c("var1", "var2", "stat", "val")
    for (row in rows) {
        for (col in cols) {
            filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_")
            if (!is.null(var1_labs)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", var1_labs[row], "_")
                var1 <- var1_labs[row]
            } else if (is.vector(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[row]), "_")
                var1 <- scientific(var1_vals[row])
            } else if (is.matrix(var1_vals)) {
                filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[1, row]), "&", scientific(var1_vals[2, row]), "_")
                var1 <- scientific(var1_vals[1, row])
            }
            if (!is.null(var2_labs)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", var2_labs[col])
                var2 <- var2_labs[col]
            } else if (is.vector(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[col]))
                var2 <- scientific(var2_vals[col])
            } else if (is.matrix(var2_vals)) {
                filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[1, col]), "&", scientific(var2_vals[2, col]))
                var2 <- scientific(var2_vals[1, col])
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



            #---Statistics: proportion of WGD
            mean_WGD_proportion <- mean(as.numeric(df_stat_sims_all$val[(which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_WGD"))]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "WGD_proportion", mean_WGD_proportion)
            #---Statistics: FGA in WGD-negative samples
            list_sims <- as.numeric(df_stat_sims_all$sim[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_WGD" & df_stat_sims_all$val == 0)])
            if (length(list_sims) == 0) {
                FGA_in_nonwgd <- NA
            } else {
                FGA_in_nonwgd <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_FGA" & df_stat_sims_all$sim %in% list_sims)]))
            }
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "FGA_in_nonwgd", FGA_in_nonwgd)
            #---Statistics: FGA in WGD-positive samples
            list_sims <- as.numeric(df_stat_sims_all$sim[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_WGD" & df_stat_sims_all$val == 1)])
            if (length(list_sims) == 0) {
                FGA_in_wgd <- NA
            } else {
                FGA_in_wgd <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_FGA" & df_stat_sims_all$sim %in% list_sims)]))
            }
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "FGA_in_wgd", FGA_in_wgd)
            #---Statistics: FGA difference between WGD-positive and WGD-negative samples
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "FGA_difference", FGA_in_wgd - FGA_in_nonwgd)
            #---Statistics: event count in WGD-negative samples
            list_sims <- as.numeric(df_stat_sims_all$sim[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_WGD" & df_stat_sims_all$val == 0)])
            if (length(list_sims) == 0) {
                event_count_in_nonwgd <- NA
            } else {
                event_count_in_nonwgd <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_event_count" & df_stat_sims_all$sim %in% list_sims)]))
            }
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "event_count_in_nonwgd", event_count_in_nonwgd)
            #---Statistics: event count in WGD-positive samples
            list_sims <- as.numeric(df_stat_sims_all$sim[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_WGD" & df_stat_sims_all$val == 1)])
            if (length(list_sims) == 0) {
                event_count_in_wgd <- NA
            } else {
                event_count_in_wgd <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "major_clone_event_count" & df_stat_sims_all$sim %in% list_sims)]))
            }
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "event_count_in_wgd", event_count_in_wgd)
            #---Statistics: event count difference between WGD-positive and WGD-negative samples
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "event_count_difference", event_count_in_wgd - event_count_in_nonwgd)



            #---Statistics: count of nonviable cells in diploid & tetraploid cells
            mean_nonviability_diploid <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=no_ploidy=2_cell_count")]))
            mean_nonviability_tetraploid <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=no_ploidy=4_cell_count")]))
            mean_nonviability_other <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=no_cell_count")])) - mean_nonviability_diploid - mean_nonviability_tetraploid
            if (is.nan(mean_nonviability_diploid)) mean_nonviability_diploid <- 0
            if (is.nan(mean_nonviability_tetraploid)) mean_nonviability_tetraploid <- 0
            if (is.nan(mean_nonviability_other)) mean_nonviability_other <- 0
            mean_nonviability_diploid <- 100 * mean_nonviability_diploid / mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=all_cell_count")]))
            mean_nonviability_tetraploid <- 100 * mean_nonviability_tetraploid / mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=all_cell_count")]))
            mean_nonviability_other <- max(0, 100 * mean_nonviability_other / mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "viability=all_cell_count")])))
            df_nonviability[nrow(df_nonviability) + 1, ] <- c(var1, var2, mean_nonviability_diploid, mean_nonviability_tetraploid, mean_nonviability_other)
            #---Statistics: count of clonal & subclonal events
            if (plot_misseg == TRUE) {
                mean_clonal_count_missegregation <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "clonal_count_missegregation")]))
                mean_subclonal_count_missegregation <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "subclonal_count_missegregation")]))
                df_missegregation[nrow(df_missegregation) + 1, ] <- c(var1, var2, mean_clonal_count_missegregation, mean_subclonal_count_missegregation)
            }
            #---Statistics: ploidy
            mean_ploidy <- mean(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy")]))
            df_stat_average[nrow(df_stat_average) + 1, ] <- c(var1, var2, "ploidy", mean_ploidy)
            #---Statistics: distribution of ploidy
            all_cell_count <- sum(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy=all_cell_count")]))
            ploidy_1_cell_count <- sum(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy=1_cell_count")]))
            ploidy_2_cell_count <- sum(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy=2_cell_count")]))
            ploidy_3_cell_count <- sum(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy=3_cell_count")]))
            ploidy_4_cell_count <- sum(as.numeric(df_stat_sims_all$val[which(df_stat_sims_all$var1 == var1 & df_stat_sims_all$var2 == var2 & df_stat_sims_all$stat == "ploidy=4_cell_count")]))
            ploidy_other_cell_count <- all_cell_count - ploidy_1_cell_count - ploidy_2_cell_count - ploidy_3_cell_count - ploidy_4_cell_count
            df_ploidy_dist[nrow(df_ploidy_dist) + 1, ] <- c(
                var1, var2,
                100 * ploidy_1_cell_count / all_cell_count,
                100 * ploidy_2_cell_count / all_cell_count,
                100 * ploidy_3_cell_count / all_cell_count,
                100 * ploidy_4_cell_count / all_cell_count,
                100 * ploidy_other_cell_count / all_cell_count
            )
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



    print(df_stat_average)



    save(df_stat_average, file = paste0(folder_workplace, "/", model_prefix, "_average_stats.rda"))
    save(df_ploidy_dist, file = paste0(folder_workplace, "/", model_prefix, "_ploidy_distribution.rda"))
    save(df_nonviability, file = paste0(folder_workplace, "/", model_prefix, "_nonviability.rda"))
    if (plot_WGD == TRUE) save(df_WGD_clonality_dist, file = paste0(folder_workplace, "/", model_prefix, "_WGD_clonality_distribution.rda"))
    if (plot_misseg == TRUE) save(df_missegregation, file = paste0(folder_workplace, "/", model_prefix, "_missegregation.rda"))
    #-----------------------------Prepare dataframe for average heatmaps
    load(paste0(folder_workplace, "/", model_prefix, "_average_stats.rda"))
    load(paste0(folder_workplace, "/", model_prefix, "_ploidy_distribution.rda"))
    load(paste0(folder_workplace, "/", model_prefix, "_nonviability.rda"))
    if (plot_WGD == TRUE) load(paste0(folder_workplace, "/", model_prefix, "_WGD_clonality_distribution.rda"))
    if (plot_misseg == TRUE) load(paste0(folder_workplace, "/", model_prefix, "_missegregation.rda"))
    df_stat_average$val <- as.numeric(df_stat_average$val)
    #   Sort dataframe according to increasing variables
    if (!is.null(var1_labs)) {
        df_stat_average$var1 <- factor(df_stat_average$var1, levels = var1_labs)
    } else if (is.vector(var1_vals)) {
        df_stat_average$var1 <- factor(scientific(as.numeric(df_stat_average$var1)), levels = scientific(sort(var1_vals)))
    } else if (is.matrix(var1_vals)) {
        df_stat_average$var1 <- factor(scientific(as.numeric(df_stat_average$var1)), levels = scientific(sort(var1_vals[1, ])))
    }
    if (!is.null(var2_labs)) {
        df_stat_average$var2 <- factor(df_stat_average$var2, levels = var2_labs)
    } else if (is.vector(var2_vals)) {
        df_stat_average$var2 <- factor(scientific(as.numeric(df_stat_average$var2)), levels = scientific(sort(var2_vals)))
    } else if (is.matrix(var2_vals)) {
        df_stat_average$var2 <- factor(scientific(as.numeric(df_stat_average$var2)), levels = scientific(sort(var2_vals[1, ])))
    }
    df_stat_average <- df_stat_average[which((is.na(df_stat_average$var1) == FALSE) & (is.na(df_stat_average$var2) == FALSE)), ]
    #--------------------------Prepare dataframe for ploidy & WGD status
    if (is.null(var1_labs)) df_ploidy_dist$var1 <- as.numeric(df_ploidy_dist$var1)
    if (is.null(var2_labs)) df_ploidy_dist$var2 <- as.numeric(df_ploidy_dist$var2)
    df_ploidy_dist$haploid <- as.numeric(df_ploidy_dist$haploid)
    df_ploidy_dist$diploid <- as.numeric(df_ploidy_dist$diploid)
    df_ploidy_dist$triploid <- as.numeric(df_ploidy_dist$triploid)
    df_ploidy_dist$tetraploid <- as.numeric(df_ploidy_dist$tetraploid)
    df_ploidy_dist$other <- as.numeric(df_ploidy_dist$other)
    #   Force identical scale for distribution maps
    if (!is.null(var1_labs)) {
        x_ticks_breaks <- var1_labs
    } else {
        x_ticks_breaks <- sort(unique(df_ploidy_dist$var1))
    }
    if (!is.null(var2_labs)) {
        y_ticks_breaks <- var2_labs
    } else {
        y_ticks_breaks <- sort(unique(df_ploidy_dist$var2))
    }
    for (ind in 1:length(x_ticks_breaks)) {
        df_ploidy_dist$var1[which(df_ploidy_dist$var1 == x_ticks_breaks[ind])] <- ind
    }
    for (ind in 1:length(y_ticks_breaks)) {
        df_ploidy_dist$var2[which(df_ploidy_dist$var2 == y_ticks_breaks[ind])] <- ind
    }
    df_ploidy_dist$var1 <- as.numeric(df_ploidy_dist$var1)
    df_ploidy_dist$var2 <- as.numeric(df_ploidy_dist$var2)
    if (plot_WGD == TRUE) {
        if (is.null(var1_labs)) df_WGD_clonality_dist$var1 <- as.numeric(df_WGD_clonality_dist$var1)
        if (is.null(var2_labs)) df_WGD_clonality_dist$var2 <- as.numeric(df_WGD_clonality_dist$var2)
        df_WGD_clonality_dist$var1_real <- df_WGD_clonality_dist$var1
        df_WGD_clonality_dist$var2_real <- df_WGD_clonality_dist$var2
        df_WGD_clonality_dist$clonal_WGD <- as.numeric(df_WGD_clonality_dist$clonal_WGD)
        df_WGD_clonality_dist$subclonal_WGD <- as.numeric(df_WGD_clonality_dist$subclonal_WGD)
        df_WGD_clonality_dist$other <- as.numeric(df_WGD_clonality_dist$other)
        #   Force identical scale for distribution maps
        for (ind in 1:length(x_ticks_breaks)) df_WGD_clonality_dist$var1[which(df_WGD_clonality_dist$var1 == x_ticks_breaks[ind])] <- ind
        for (ind in 1:length(y_ticks_breaks)) df_WGD_clonality_dist$var2[which(df_WGD_clonality_dist$var2 == y_ticks_breaks[ind])] <- ind
        df_WGD_clonality_dist$var1 <- as.numeric(df_WGD_clonality_dist$var1)
        df_WGD_clonality_dist$var2 <- as.numeric(df_WGD_clonality_dist$var2)
    }
    #-------------------------Prepare dataframe for nonviable cell count
    if (is.null(var1_labs)) df_nonviability$var1 <- as.numeric(df_nonviability$var1)
    if (is.null(var2_labs)) df_nonviability$var2 <- as.numeric(df_nonviability$var2)
    df_nonviability$var1_real <- df_nonviability$var1
    df_nonviability$var2_real <- df_nonviability$var2
    df_nonviability$diploid <- as.numeric(df_nonviability$diploid)
    df_nonviability$tetraploid <- as.numeric(df_nonviability$tetraploid)
    df_nonviability$other <- as.numeric(df_nonviability$other)
    df_nonviability$total <- df_nonviability$diploid + df_nonviability$tetraploid + df_nonviability$other
    #   Force identical scale for distribution maps
    for (ind in 1:length(x_ticks_breaks)) df_nonviability$var1[which(df_nonviability$var1 == x_ticks_breaks[ind])] <- ind
    for (ind in 1:length(y_ticks_breaks)) df_nonviability$var2[which(df_nonviability$var2 == y_ticks_breaks[ind])] <- ind
    df_nonviability$var1 <- as.numeric(df_nonviability$var1)
    df_nonviability$var2 <- as.numeric(df_nonviability$var2)
    #---------------------------------Prepare dataframe for event counts
    if (plot_misseg == TRUE) {
        if (is.null(var1_labs)) df_missegregation$var1 <- as.numeric(df_missegregation$var1)
        if (is.null(var2_labs)) df_missegregation$var2 <- as.numeric(df_missegregation$var2)
        df_missegregation$var1_real <- df_missegregation$var1
        df_missegregation$var2_real <- df_missegregation$var2
        df_missegregation$clonal_count_missegregation <- as.numeric(df_missegregation$clonal_count_missegregation)
        df_missegregation$subclonal_count_missegregation <- as.numeric(df_missegregation$subclonal_count_missegregation)
        df_missegregation$total_count_missegregation <- df_missegregation$clonal_count_missegregation + df_missegregation$subclonal_count_missegregation
        df_missegregation$Clonal <- 100 * df_missegregation$clonal_count_missegregation / df_missegregation$total_count_missegregation
        df_missegregation$Subclonal <- 100 * df_missegregation$subclonal_count_missegregation / df_missegregation$total_count_missegregation
        #   Force identical scale for distribution maps
        for (ind in 1:length(x_ticks_breaks)) df_missegregation$var1[which(df_missegregation$var1 == x_ticks_breaks[ind])] <- ind
        for (ind in 1:length(y_ticks_breaks)) df_missegregation$var2[which(df_missegregation$var2 == y_ticks_breaks[ind])] <- ind
        df_missegregation$var1 <- as.numeric(df_missegregation$var1)
        df_missegregation$var2 <- as.numeric(df_missegregation$var2)
    }
    if (is.null(var1_labs)) {
        x_ticks_breaks_val <- x_ticks_breaks
        x_ticks_breaks <- scientific(x_ticks_breaks)
    } else {
        x_ticks_breaks_val <- var1_labs
    }
    if (is.null(var2_labs)) {
        y_ticks_breaks_val <- y_ticks_breaks
        if (is.null(var2_labs)) y_ticks_breaks <- scientific(y_ticks_breaks)
    } else {
        y_ticks_breaks_val <- var2_labs
    }
    #---------------------------------------Plot statistics: clone count
    filename <- paste0(model_prefix, "_1_clone_count.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "clone_count"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "RdPu", name = "Clone count") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #---------------------------Plot statistics: Shannon diversity index
    filename <- paste0(model_prefix, "_2_Shannon_index.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "Shannon_index"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_viridis(discrete = FALSE, name = "Shannon index") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #-------------Plot statistics: fitness with respect to diploid clone
    df_stat_average$val[which(df_stat_average$stat == "fitness")] <- log(df_stat_average$val[which(df_stat_average$stat == "fitness")])
    filename <- paste0(model_prefix, "_3_fitness.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "fitness"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "YlOrBr", name = "Log(fitness)") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #--------------------------------------------Plot statistics: ploidy
    filename <- paste0(model_prefix, "_4_ploidy.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "ploidy"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_gradientn(
            colours = c("azure4", "royalblue3", "ghostwhite", "seagreen2", "tomato1", "azure4"),
            breaks = c(0, 1, 2, 3, 4, 5),
            labels = c(0, 1, 2, 3, 4, 5),
            limits = c(0, 5),
            name = "Average ploidy"
        ) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
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
        labs(fill = "Ploidy") +
        xlab(var1_lab) +
        scale_x_continuous(breaks = 1:length(x_ticks_breaks), labels = x_ticks_breaks) +
        ylab(var2_lab) +
        scale_y_continuous(breaks = 1:length(y_ticks_breaks), labels = y_ticks_breaks) +
        scale_fill_manual(
            values = c(
                "haploid" = "royalblue3",
                "diploid" = "ghostwhite",
                "triploid" = "seagreen2",
                "tetraploid" = "tomato1",
                "other" = "azure4"
            ),
            breaks = c("haploid", "diploid", "triploid", "tetraploid", "other"),
            labels = c("1", "2", "3", "4", ">4")
        ) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2, "cm"))
    print(p)
    dev.off()
    #----------------Plot statistics: count of clonal & subclonal events
    custom_geom_scatterpie_legend <- function(radius, x, y, n = 5, labeller, textsize = 1) {
        if (length(radius) > n) {
            radius <- unique(sapply(seq(min(radius), max(radius),
                length.out = n
            ), scatterpie:::round_digit))
        }
        label <- FALSE
        if (!missing(labeller)) {
            if (!inherits(labeller, "function")) {
                stop("labeller should be a function for converting radius")
            }
            label <- TRUE
        }
        dd <- data.frame(
            r = radius, start = 0, end = 2 * pi, x = x,
            y = y + radius - max(radius), maxr = max(radius)
        )
        if (label) {
            dd$label <- labeller(dd$r)
        } else {
            dd$label <- dd$r
        }
        list(
            ggforce:::geom_arc_bar(aes_(
                x0 = ~x, y0 = ~y, r0 = ~r, r = ~r,
                start = ~start, end = ~end
            ), data = dd, inherit.aes = FALSE),
            geom_segment(aes_(x = ~x, xend = ~ x + maxr * 1.5, y = ~ y +
                r, yend = ~ y + r), data = dd, inherit.aes = FALSE),
            geom_text(aes_(x = ~ x + maxr * 1.6, y = ~ y + r, label = ~label),
                data = dd, hjust = "left", inherit.aes = FALSE, size = textsize
            )
        )
    }
    if (plot_misseg == TRUE) {
        filename <- paste0(model_prefix, "_5_count_missegregations.jpeg")
        jpeg(file = filename, width = 1000, height = 1100)
        max_total_count_missegregation <- max(df_missegregation$total_count_missegregation)
        p <- ggplot() +
            geom_scatterpie(aes(x = var1, y = var2, r = 0.45 * total_count_missegregation / max_total_count_missegregation),
                data = df_missegregation,
                cols = c("Clonal", "Subclonal")
            ) +
            coord_equal() +
            labs(fill = "Missegregations") +
            xlab(var1_lab) +
            scale_x_continuous(limits = c(0, 1 + max(df_missegregation$var1)), expand = c(0, 0), breaks = 1:length(x_ticks_breaks), labels = x_ticks_breaks) +
            ylab(var2_lab) +
            scale_y_continuous(limits = c(0, 1 + max(df_missegregation$var2)), expand = c(0, 0), breaks = 1:length(y_ticks_breaks), labels = y_ticks_breaks) +
            scale_fill_manual(values = c("orchid4", "moccasin")) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
        for (i in seq(1, length(y_ticks_breaks), by = 1)) {
            var1_all <- df_missegregation$var1_real
            if (is.null(var1_labs)) var1_all <- scientific(df_missegregation$var1_real)
            var2_all <- df_missegregation$var2_real
            if (is.null(var2_labs)) var2_all <- scientific(df_missegregation$var2_real)
            tmp <- df_missegregation$total_count_missegregation[which(df_missegregation$var1_real == x_ticks_breaks_val[length(x_ticks_breaks_val)] & df_missegregation$var2_real == y_ticks_breaks_val[i])]
            if (tmp <= 0) next
            p <- p + custom_geom_scatterpie_legend(0.45 * tmp / max_total_count_missegregation, x = length(x_ticks_breaks_val), y = i, labeller = function(x) round(x * max_total_count_missegregation / 0.45), textsize = 7)
        }
        print(p)
        dev.off()
    }
    #----------------Plot statistics: distribution of WGD (sub)clonality
    if (plot_WGD == TRUE) {
        filename <- paste0(model_prefix, "_6_WGD_clonality_distribution.jpeg")
        jpeg(file = filename, width = 1000, height = 1100)
        p <- ggplot() +
            geom_scatterpie(aes(x = var1, y = var2, r = 0.45),
                data = df_WGD_clonality_dist,
                cols = c("clonal_WGD", "subclonal_WGD", "other")
            ) +
            coord_equal() +
            labs(fill = "WGD") +
            xlab(var1_lab) +
            scale_x_continuous(breaks = 1:length(x_ticks_breaks), labels = x_ticks_breaks) +
            ylab(var2_lab) +
            scale_y_continuous(breaks = 1:length(y_ticks_breaks), labels = y_ticks_breaks) +
            scale_fill_manual(
                values = c(
                    "clonal_WGD" = "sienna",
                    "subclonal_WGD" = "wheat",
                    "other" = "gray"
                ),
                breaks = c("clonal_WGD", "subclonal_WGD", "other"),
                labels = c("Clonal", "Subclonal", "Other")
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
        print(p)
        dev.off()
    }
    #---------------------Plot statistics: event count before clonal WGD
    if (plot_WGD == TRUE) {
        filename <- paste0(model_prefix, "_7_event_count_pre_clonal_WGD.jpeg")
        jpeg(file = filename, width = 1000, height = 1100)
        df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "event_count_pre_clonal_WGD"), ]
        p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
            geom_tile() +
            scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
            scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
            coord_equal() +
            xlab(var1_lab) +
            ylab(var2_lab) +
            scale_fill_gradientn(
                colours = c("white", "chartreuse4"),
                na.value = "black",
                name = "Event count before clonal WGD"
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(1.5, "cm"))
        print(p)
        dev.off()
    }
    #-----------------Plot statistics: nonviability percentage by ploidy
    filename <- paste0(model_prefix, "_8_count_nonviability.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    max_total <- max(df_nonviability$total)
    p <- ggplot() +
        geom_scatterpie(aes(x = var1, y = var2, r = 0.45 * total / max_total),
            data = df_nonviability,
            cols = c("diploid", "tetraploid", "other")
        ) +
        coord_equal() +
        labs(fill = "Nonviability") +
        xlab(var1_lab) +
        scale_x_continuous(limits = c(0, 1 + max(df_nonviability$var1)), expand = c(0, 0), breaks = 1:length(x_ticks_breaks), labels = x_ticks_breaks) +
        ylab(var2_lab) +
        scale_y_continuous(limits = c(0, 1 + max(df_nonviability$var2)), expand = c(0, 0), breaks = 1:length(y_ticks_breaks), labels = y_ticks_breaks) +
        # scale_fill_manual(values = c("orchid4", "moccasin")) +
        scale_fill_manual(
            values = c(
                "diploid" = "ghostwhite",
                "tetraploid" = "tomato1",
                "other" = "azure4"
            ),
            breaks = c("diploid", "tetraploid", "other"),
            labels = c("Diploid", "Tetraploid", "Other")
        ) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    for (i in seq(1, length(y_ticks_breaks), by = 1)) {
        var1_all <- df_nonviability$var1_real
        if (is.null(var1_labs)) var1_all <- scientific(df_nonviability$var1_real)
        var2_all <- df_nonviability$var2_real
        if (is.null(var2_labs)) var2_all <- scientific(df_nonviability$var2_real)
        tmp <- df_nonviability$total[which(df_nonviability$var1_real == x_ticks_breaks_val[length(x_ticks_breaks_val)] & df_nonviability$var2_real == y_ticks_breaks_val[i])]
        if (tmp <= 0) next
        p <- p + custom_geom_scatterpie_legend(0.45 * tmp / max_total, x = length(x_ticks_breaks_val), y = i, labeller = function(x) paste0(round(x * max_total / 0.45), "%"), textsize = 7)
    }
    print(p)
    dev.off()




    #-----------Plot statistics: Proportion of WGD-dominated simulations
    filename <- paste0(model_prefix, "_9_WGD_proportion.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "WGD_proportion"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "YlOrBr", name = "WGD proportion") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #--Plot statistics: FGA difference between WGD & non-WGD simulations
    filename <- paste0(model_prefix, "_10_WGD_FGA_difference.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "FGA_difference"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "YlOrBr", name = "WGD FGA difference") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
    #--Plot statistics: event count difference between WGD & non-WGD simulations
    filename <- paste0(model_prefix, "_11_WGD_aneuploidy_difference.jpeg")
    jpeg(file = filename, width = 1000, height = 1100)
    df_stat_average_plot <- df_stat_average[which(df_stat_average$stat == "event_count_difference"), ]
    p <- ggplot(df_stat_average_plot, aes(var1, var2, fill = val)) +
        geom_tile() +
        scale_x_discrete(expand = c(1 / length(rows), 1 / length(rows))) +
        scale_y_discrete(expand = c(1 / length(cols), 1 / length(cols))) +
        coord_equal() +
        xlab(var1_lab) +
        ylab(var2_lab) +
        scale_fill_distiller(palette = "YlOrBr", name = "WGD aneuploidy difference") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 40), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()
}

#' @export
statistics_multivar_one_simulation <- function(filename, var1_name, var2_name, var1, var2, sim, plot_WGD, plot_misseg) {
    load(filename)
    #--------------------------------Create dataframe for all statistics
    df_stat_sim <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_stat_sim) <- c("var1", "var2", "sim", "stat", "val")
    #---------------------------Input required variables from simulation
    genotype_list_selection_rate <- simulation$clonal_evolution$genotype_list_selection_rate
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    genotype_list_WGD_count <- simulation$clonal_evolution$genotype_list_WGD_count
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
    N_chromosomes <<- length(genotype_list_ploidy_chrom[[1]])
    vec_CN_block_no <<- rep(0, N_chromosomes)
    for (chrom in 1:N_chromosomes) {
        vec_CN_block_no[chrom] <<- length(genotype_list_ploidy_block[[1]][[chrom]][[1]])
    }
    vec_chromosome_id <<- 1:N_chromosomes
    size_CN_block_DNA <<- 500000
    #-------------------Function to find clonal ancestry of given clones
    find_clonal_ancestry <- function(list_subclonal_ancestry) {
        if (length(list_subclonal_ancestry) == 0) {
            clonal_ancestry <- c()
        } else if (length(list_subclonal_ancestry) == 1) {
            clonal_ancestry <- list_subclonal_ancestry[[1]]
        } else {
            clonal_ancestry <- list_subclonal_ancestry[[1]]
            for (i in 2:length(list_subclonal_ancestry)) {
                clonal_ancestry <- intersect(clonal_ancestry, list_subclonal_ancestry[[i]])
            }
        }
        return(clonal_ancestry)
    }
    #-------------------------------------------------Find unique clones
    table_clone <- as.data.frame(table(ID = simulation$sample$sample_clone_ID))
    table_clone$ID <- as.numeric(as.vector(table_clone$ID))
    table_clone$ID_unique <- 0
    table_clone$ID_unique[1] <- 1
    ID_unique <- 1
    if (nrow(table_clone) > 1) {
        for (j in 2:nrow(table_clone)) {
            clone_new <- table_clone$ID[j]
            ID_unique_new <- 0
            for (i in 1:(j - 1)) {
                clone_old <- table_clone$ID[i]
                #   Check if selection rate is same
                if (genotype_list_selection_rate[clone_new] != genotype_list_selection_rate[clone_old]) next
                #   Check if CN profile is same
                if (any(genotype_list_ploidy_chrom[[clone_new]] != genotype_list_ploidy_chrom[[clone_old]])) next
                tmp <- 1
                for (chrom in 1:length(genotype_list_ploidy_chrom[[clone_new]])) {
                    if (genotype_list_ploidy_chrom[[clone_new]][chrom] <= 0) next
                    for (strand in 1:genotype_list_ploidy_chrom[[clone_new]][chrom]) {
                        if (!setequal(genotype_list_ploidy_block[[clone_new]][[chrom]][[strand]], genotype_list_ploidy_block[[clone_old]][[chrom]][[strand]])) tmp <- 0
                    }
                }
                if (tmp == 0) next
                ID_unique_new <- table_clone$ID_unique[i]
            }
            if (ID_unique_new == 0) {
                ID_unique <- ID_unique + 1
                ID_unique_new <- ID_unique
            }
            table_clone$ID_unique[j] <- ID_unique_new
        }
    }
    table_clone_unique <- data.frame(ID_unique = 1:ID_unique)
    table_clone_unique$Freq <- 0
    for (ID_unique in 1:nrow(table_clone_unique)) {
        table_clone_unique$Freq[ID_unique] <- sum(table_clone$Freq[which(table_clone$ID_unique == ID_unique)])
    }
    Clone_ID <- table_clone$ID
    Clone_ID_max <- Clone_ID[which(table_clone$Freq == max(table_clone$Freq))][1]
    Clone_ID_unique <- as.numeric(as.vector(table_clone_unique$ID_unique))
    #---------------------------------------Find ancestry of every clone
    subclonal_ancestry <- vector("list", length(Clone_ID))
    for (i in 1:length(Clone_ID)) {
        ancestry <- Clone_ID[i]
        while (ancestry[1] != 0) ancestry <- c(evolution_origin[ancestry[1]], ancestry)
        subclonal_ancestry[[i]] <- ancestry
    }
    #-------------------------------------------------Statistics: ploidy
    compute_ploidy <- function(vec_CN_block_no, ploidy_chrom, ploidy_block) {
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
    #   Find ploidy and selection rate for each unique clone and its mother clone
    table_clone$ploidy <- 0
    table_clone$mother_ploidy <- 0
    for (i in 1:nrow(table_clone)) {
        table_clone$ploidy[i] <- compute_ploidy(vec_CN_block_no, genotype_list_ploidy_chrom[[Clone_ID[i]]], genotype_list_ploidy_block[[Clone_ID[i]]])
        if (evolution_origin[Clone_ID[i]] <= 0) {
            mother_clone <- Clone_ID[i]
        } else {
            mother_clone <- evolution_origin[Clone_ID[i]]
        }
        table_clone$mother_ploidy[i] <- compute_ploidy(vec_CN_block_no, genotype_list_ploidy_chrom[[mother_clone]], genotype_list_ploidy_block[[mother_clone]])
    }
    table_clone$rounded_ploidy <- round(table_clone$ploidy)
    table_clone$rounded_mother_ploidy <- round(table_clone$mother_ploidy)
    ploidy_unique <- unique(table_clone$rounded_ploidy)
    #   Compute mean ploidy for simulation
    ploidy <- sum((table_clone$ploidy) * (table_clone$Freq)) / sum(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "ploidy", ploidy)
    #-----------------------------Statistics: cell count for each ploidy
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "ploidy=all_cell_count", sum(table_clone$Freq))
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_cell_count"), sum(table_clone$Freq[which(table_clone$rounded_ploidy == ploidy)]))
    }
    #--------------------------------------------Statistics: clone count
    clone_count <- length(Clone_ID_unique)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clone_count", clone_count)
    #------------------Statistics: fitness with respect to diploid clone
    table_clone$fitness <- genotype_list_selection_rate[Clone_ID] / genotype_list_selection_rate[1]
    fitness <- sum((table_clone$fitness) * (table_clone$Freq)) / sum(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "fitness", fitness)
    #-------------------Statistics: cell count for each viability status
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viability=all_cell_count", sum(table_clone$Freq))
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viability=yes_cell_count", sum(table_clone$Freq[which(table_clone$fitness > 0)]))
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viability=no_cell_count", sum(table_clone$Freq[which(table_clone$fitness <= 0)]))
    #------------Statistics: cell count for nonviability for each ploidy
    mother_ploidy_unique <- c(2, 4)
    for (i in 1:length(mother_ploidy_unique)) {
        ploidy <- mother_ploidy_unique[i]
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viability=no_ploidy=", ploidy, "_cell_count"), sum(table_clone$Freq[which((table_clone$fitness <= 0) & (table_clone$rounded_mother_ploidy == ploidy))]))
    }
    #--------------------------------Statistics: Shannon diversity index
    Shannon_index <- diversity(table_clone_unique$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "Shannon_index", Shannon_index)
    #--------------------------------Statistics: main clone's WGD status
    if (genotype_list_WGD_count[Clone_ID_max] > 0) {
        Clone_ID_WGD_status <- 1
    } else {
        Clone_ID_WGD_status <- 0
    }
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "major_clone_WGD", Clone_ID_WGD_status)
    #----------------Statistics: main clone's Fraction of Genome Altered
    get_FGA <- function(package_clonal_evolution,
                        clone_ID) {
        ploidy_normalization <- TRUE
        plotcol <- "state"
        fillna <- TRUE

        WGD_count <- package_clonal_evolution$genotype_list_WGD_count[clone_ID]
        state_mode <- 2^(WGD_count + 1)

        #   Find total CN profile
        CNbins_sims <- get_cn_profile(package_clonal_evolution, clone_ID)
        CNbins_sims$cell_id <- paste0("SIMULATION1-Library-1-1")
        CNbins_sims$chr <- as.character(CNbins_sims$chr)
        class(CNbins_sims) <- "data.frame"
        copynumber_sims <- createCNmatrix(CNbins_sims,
            field = plotcol, wholegenome = FALSE,
            fillnaplot = fillna, centromere = FALSE
        )
        if (ploidy_normalization == TRUE) {
            copynumber_sims <- normalize_cell_ploidy(copynumber_sims, state_mode, round = FALSE)
        }
        #   Statistics - FGA
        sample_CN <- copynumber_sims[[paste0("SIMULATION1-Library-1-1")]]
        FGA <- length(which(sample_CN != 2)) / length(sample_CN)
        #####
        #####
        #####
        #####
        #####
        # if (WGD_count > 0) {
        #     genotype_list_ploidy_chrom <- package_clonal_evolution$genotype_list_ploidy_chrom
        #     evolution_genotype_changes <- package_clonal_evolution$evolution_genotype_changes
        #     evolution_traj_time <- package_clonal_evolution$evolution_traj_time
        #     evolution_traj_clonal_ID <- package_clonal_evolution$evolution_traj_clonal_ID
        #     evolution_origin <- package_clonal_evolution$evolution_origin
        #     pre_misseg <- 0
        #     pre_arm_misseg <- 0
        #     WGD_clones <- c()
        #     post_misseg <- 0
        #     post_arm_misseg <- 0
        #     clone <- clone_ID
        #     while (clone > 0) {
        #         n_misseg <- 0
        #         n_arm_misseg <- 0
        #         genotype_changes <- evolution_genotype_changes[[clone]]
        #         if (length(genotype_changes) > 0) {
        #             for (i in 1:length(genotype_changes)) {
        #                 if (genotype_changes[[i]][1] == "missegregation") {
        #                     n_misseg <- n_misseg + 1
        #                 } else if (genotype_changes[[i]][1] == "chromosome-arm-missegregation") {
        #                     post_arm_misseg <- post_arm_misseg + 1
        #                 } else if (genotype_changes[[i]][1] == "whole-genome-duplication") {
        #                     WGD_clones <- c(WGD_clones, clone)
        #                 }
        #             }
        #         }
        #         if (length(WGD_clones) > 0) {
        #             post_misseg <- post_misseg + n_misseg
        #             post_arm_misseg <- post_arm_misseg + n_arm_misseg
        #         } else {
        #             pre_misseg <- pre_misseg + n_misseg
        #             pre_arm_misseg <- pre_arm_misseg + n_arm_misseg
        #         }
        #         clone <- evolution_origin[clone]
        #     }
        #     for (i in 1:length(evolution_traj_time)) {
        #         if (WGD_clones %in% evolution_traj_clonal_ID[[i]]) {
        #             WGD_age <- evolution_traj_time[i]
        #             break
        #         }
        #     }
        #     cat("-------------------------------------------------------\n")
        #     cat(paste0("Pre-WGD missegregations      = ", pre_misseg, "\n"))
        #     cat(paste0("Pre-WGD arm-missegregations  = ", pre_arm_misseg, "\n"))
        #     cat(paste0("WGD ancestor                 = ", WGD_clones, "\n"))
        #     cat(paste0("WGD ancestor age             = ", WGD_age / 365, "\n"))
        #     cat(paste0("WGD ancestor CN profile      : \n"))
        #     print(genotype_list_ploidy_chrom[[WGD_clones]])
        #     cat(paste0("Post-WGD missegregations     = ", post_misseg, "\n"))
        #     cat(paste0("Post-WGD arm-missegregations = ", post_arm_misseg, "\n"))
        #     cat(paste0("Final CN profile             : \n"))
        #     print(genotype_list_ploidy_chrom[[clone_ID]])
        #     cat("-------------------------------------------------------\n")
        # }
        #####
        #####
        #####
        #####
        #####
        #   Output statistics
        return(FGA)
    }
    Clone_ID_FGA <- get_FGA(simulation$clonal_evolution, Clone_ID_max)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "major_clone_FGA", Clone_ID_FGA)
    #-------------------------------Statistics: main clone's event count
    get_event_count <- function(package_clonal_evolution,
                                clone_ID) {
        event_count <- 0
        clone <- clone_ID
        while (clone > 0) {
            genotype_changes <- evolution_genotype_changes[[clone]]
            event_count <- event_count + length(genotype_changes)
            clone <- evolution_origin[clone]
        }
        return(event_count)
    }
    Clone_ID_event_count <- get_event_count(simulation$clonal_evolution, Clone_ID_max)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "major_clone_event_count", Clone_ID_event_count)
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
    #-----------------------Statistics: count of clonal/subclonal events
    #-----------------------------------------------for all alive clones
    find_event_count <- function(ancestry, event_type) {
        event_count <- 0
        if (length(ancestry) > 0) {
            for (i in 1:length(ancestry)) {
                if (ancestry[i] == 0) next
                if (length(evolution_genotype_changes[[ancestry[i]]]) == 0) next
                for (j in 1:length(evolution_genotype_changes[[ancestry[i]]])) {
                    if (evolution_genotype_changes[[ancestry[i]]][[j]][1] == event_type) {
                        event_count <- event_count + 1
                    }
                }
            }
        }
        return(event_count)
    }
    if (plot_misseg == TRUE) {
        #   Find count of clonal missegs
        clonal_ancestry <- find_clonal_ancestry(subclonal_ancestry)
        clonal_count_missegregation <- find_event_count(clonal_ancestry, "missegregation")
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clonal_count_missegregation", clonal_count_missegregation)
        #   Find average count of subclonal missegs
        table_clone$count_missegregation <- 0
        for (k in 1:nrow(table_clone)) table_clone$count_missegregation[k] <- find_event_count(subclonal_ancestry[[k]], "missegregation")
        subclonal_count_missegregation <- sum((table_clone$count_missegregation) * (table_clone$Freq)) / sum(table_clone$Freq) - clonal_count_missegregation
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "subclonal_count_missegregation", subclonal_count_missegregation)
    }
    #-----------------------Statistics: count of clonal/subclonal events
    #-----------------------------------------for clones based on ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        mini_table_clone <- table_clone[which(table_clone$rounded_ploidy == ploidy), ]
        mini_subclonal_ancestry <- subclonal_ancestry[which(table_clone$rounded_ploidy == ploidy)]
        #   Find count of clonal missegs
        clonal_ancestry <- find_clonal_ancestry(mini_subclonal_ancestry)
        clonal_count_missegregation <- find_event_count(clonal_ancestry, "missegregation")
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_clonal_count_missegregation"), clonal_count_missegregation)
        #   Find average count of subclonal missegs
        subclonal_count_missegregation <- sum((mini_table_clone$count_missegregation) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_subclonal_count_missegregation"), subclonal_count_missegregation)
    }
    #--------------------------Statistics: event count before clonal WGD
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

#' @export
statistics_multivar_vector <- function(model_prefix = "",
                                       folder_workplace = NULL,
                                       folder_plots = NULL,
                                       model_variables_base,
                                       var1_name = "",
                                       var1_vals = c(),
                                       var1_labs = NULL,
                                       var2_name = "",
                                       var2_vals = c(),
                                       var2_labs = NULL,
                                       var_labs = NULL,
                                       name_lab = "",
                                       plotname = "",
                                       plot_WGD = FALSE,
                                       n_simulations = 0,
                                       compute_parallel = TRUE,
                                       example_simulation = FALSE,
                                       n_cores = NULL,
                                       R_libPaths = NULL) {
    library(scales)
    library(data.table)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    if (!is.null(var1_labs)) {
        inds <- 1:length(var1_labs)
    } else if (is.vector(var1_vals)) {
        inds <- 1:length(var1_vals)
    } else if (is.matrix(var1_vals)) {
        inds <- 1:ncol(var1_vals)
    }
    #------------------------------------Get statistics from simulations
    df_stat_sims_all_list <- vector("list", length(inds))
    for (ind in inds) {
        cat(paste("\nSTATISTICS FOR BATCH ", ind, "/", length(inds), "...\n", sep = ""))
        filename_prefix <<- paste0(folder_workplace, "/", model_prefix, "_")
        if (!is.null(var1_labs)) {
            filename_prefix <<- paste0(filename_prefix, var1_name, "=", var1_labs[ind], "_")
            var1 <- var1_labs[ind]
        } else if (is.vector(var1_vals)) {
            filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[ind]), "_")
            var1 <- scientific(var1_vals[ind])
        } else if (is.matrix(var1_vals)) {
            filename_prefix <<- paste0(filename_prefix, var1_name, "=", scientific(var1_vals[1, ind]), "&", scientific(var1_vals[2, ind]), "_")
            var1 <- scientific(var1_vals[1, ind])
        }
        if (!is.null(var2_labs)) {
            filename_prefix <<- paste0(filename_prefix, var2_name, "=", var2_labs[ind])
            var2 <- var2_labs[ind]
        } else if (is.vector(var2_vals)) {
            filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[ind]))
            var2 <- scientific(var2_vals[ind])
        } else if (is.matrix(var2_vals)) {
            filename_prefix <<- paste0(filename_prefix, var2_name, "=", scientific(var2_vals[1, ind]), "&", scientific(var2_vals[2, ind]))
            var2 <- scientific(var2_vals[1, ind])
        }
        if (compute_parallel == FALSE) {
            #--Get statistics for each simulation in sequential mode
            df_stat_sims_list <- vector("list", n_simulations)
            for (sim in 1:n_simulations) {
                filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
                df_stat_sim <- simulation_statistics_one_simulation(filename, var1_name, var2_name, var1, var2, sim, plot_WGD)
                df_stat_sims_list[[sim]] <- df_stat_sim
            }
        } else {
            #----Get statistics for each simulation in parallel mode
            library(parallel)
            library(pbapply)
            start_time <- Sys.time()
            #   Start parallel cluster
            if (is.null(n_cores)) {
                numCores <- detectCores()
            } else {
                numCores <- n_cores
            }
            cl <- makePSOCKcluster(numCores - 1)
            if (is.null(R_libPaths) == FALSE) {
                R_libPaths <<- R_libPaths
                clusterExport(cl, varlist = c(
                    "R_libPaths"
                ))
                clusterEvalQ(cl = cl, .libPaths(R_libPaths))
            }
            clusterEvalQ(cl = cl, library(scales))
            clusterEvalQ(cl = cl, library(vegan))
            #   Prepare input parameters for plotting
            var1_name <<- var1_name
            var2_name <<- var2_name
            var1 <<- var1
            var2 <<- var2
            plot_WGD <<- plot_WGD
            simulation_statistics_one_simulation <<- simulation_statistics_one_simulation
            clusterExport(cl, varlist = c(
                "filename_prefix",
                "var1_name",
                "var2_name",
                "var1",
                "var2",
                "plot_WGD",
                "simulation_statistics_one_simulation"
            ))
            #   Get statistics in parallel
            pbo <- pboptions(type = "txt")
            df_stat_sims_list <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(sim) {
                filename <- paste0(filename_prefix, "_simulation_", sim, ".rda")
                df_stat_sim <- simulation_statistics_one_simulation(filename, var1_name, var2_name, var1, var2, sim, plot_WGD)
                return(df_stat_sim)
            })
            #   Stop parallel cluster
            stopCluster(cl)
            end_time <- Sys.time()
            print(end_time - start_time)
        }
        df_tmp <- rbindlist(df_stat_sims_list)
        df_tmp$val <- as.numeric(df_tmp$val)
        df_tmp$sim <- as.numeric(df_tmp$sim)
        df_tmp$lab <- var_labs[ind]
        df_stat_sims_all_list[[ind]] <- df_tmp
        #--------------------------Make plots for individual simulations
        if (example_simulation == TRUE) {
            #   Choose simulation closest to median clonal CCF score
            vec_sim_score <- rep(0, length(unique(df_plot$sim)))
            for (sim in 1:length(vec_sim_score)) {
                for (group in 1:5) {
                    vec_sim_score[sim] <- vec_sim_score[sim] + 10^(3 * (5 - group)) * df_plot$val[which(df_plot$sim == sim & df_plot$stat == as.character(group))]
                }
            }
            loc <- which.min(abs(vec_sim_score - median(vec_sim_score)))[1]
            #   Plot CN profile for this simulation
            plot_cn_heatmap_one_simulation(
                model = plotname_prefix,
                loc,
                folder_workplace = paste0(folder_workplace, "/"),
                folder_plots = "",
                plotcol = "total-copy",
                plottree = FALSE,
                plotfrequency = FALSE,
                show_legend = FALSE,
                show_library_label = FALSE,
                show_clone_label = FALSE,
                CN_data = "TRUTH",
                phylo = "TRUTH",
                filename_suffix = "_EXAMPLE",
                width = 1000,
                height = 1000
            )
        }
    }
    df_stat_sims_all <- rbindlist(df_stat_sims_all_list)
    save(df_stat_sims_all, file = paste0(folder_workplace, "/", plotname, ".rda"))
    #---------------------------------------------Functions for plotting
    load(paste0(folder_workplace, "/", plotname, ".rda"))
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
    geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                                  draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                                  show.legend = NA, inherit.aes = TRUE) {
        layer(
            data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
            position = position, show.legend = show.legend, inherit.aes = inherit.aes,
            params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
        )
    }
    #   Function to plot MRCA age
    plot_MRCA_age <- function(df_plot, lab) {
        p_MRCA <- ggplot(df_plot, aes(x = lab, y = -val)) +
            geom_violin(fill = "azure4", width = 1, alpha = 0.2) +
            geom_boxplot(fill = "azure4", width = 0.1) +
            coord_flip() +
            xlab(name_lab) +
            ylab(lab) +
            scale_y_continuous(breaks = c(-1, -0.5, 0), labels = c("-1", "-0.5", "0"), limits = c(-1, 0)) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = "top")
        return(p_MRCA)
    }
    #   Function to plot clonal or subclonal events
    plot_events <- function(df_plot, lab, pos_legend, limits = NULL) {
        p_events <- ggplot(df_plot, aes(x = lab, y = val, fill = stat)) +
            geom_split_violin(width = 1, alpha = 0.2) +
            geom_boxplot(width = 0.1) +
            coord_flip() +
            labs(fill = "") +
            ylab(lab) +
            scale_fill_manual(
                values = c(
                    "gain" = "indianred3",
                    "loss" = "dodgerblue3"
                ),
                breaks = c("gain", "loss"),
                labels = c("Gain", "Loss")
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = pos_legend, axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 0))
        if (!is.null(limits)) {
            p_events <- p_events + ylim(limits)
        }
        return(p_events)
    }
    #   Function to plot both clonal and subclonal events
    plot_clonal_and_subclonal_events <- function(df_plot, lab, var_labs, pos_legend, limits = NULL) {
        p_events <- ggplot(df_plot, aes(x = lab, y = val, fill = stat)) +
            geom_boxplot(width = 0.4) +
            coord_flip() +
            labs(fill = "") +
            ylab(lab) +
            scale_fill_manual(
                values = c(
                    "clonal_gain" = "#ff0000",
                    "clonal_loss" = "#0072B2",
                    "subclonal_gain" = "olivedrab",
                    "subclonal_loss" = "gold2"
                ),
                breaks = c("clonal_gain", "clonal_loss", "subclonal_gain", "subclonal_loss"),
                labels = c("Clonal gain", "Clonal loss", "Subclonal gain", "Subclonal loss")
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = pos_legend, legend.title = element_text(size = 0), legend.text = element_text(size = 20), axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 0)) +
            guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        if (!is.null(limits)) {
            p_events <- p_events + ylim(limits)
        }
        return(p_events)
    }
    #   Function to plot fitness for both diploid and tetraploid cells
    plot_fitness_by_ploidy <- function(df_plot, lab, var_labs, pos_legend, limits = NULL) {
        p_events <- ggplot(df_plot, aes(x = lab, y = val, fill = stat)) +
            # geom_violin(width = 1, alpha = 0.2) +
            geom_split_violin(width = 1.5, alpha = 0.2) +
            geom_boxplot(width = 0.2) +
            coord_flip() +
            labs(fill = "") +
            ylab(lab) +
            scale_fill_manual(
                values = c(
                    "ploidy=2" = "magenta4",
                    "ploidy=4" = "darkorange"
                ),
                breaks = c("ploidy=2", "ploidy=4"),
                labels = c("ploidy=2", "ploidy=4")
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = pos_legend, legend.title = element_text(size = 0), legend.text = element_text(size = 20), axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 0)) +
            guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        if (!is.null(limits)) {
            p_events <- p_events + ylim(limits)
        }
        return(p_events)
    }
    #   Function to plot fitness for both WGD and non-WGD cells
    plot_fitness_by_WGD_status <- function(df_plot, lab, var_labs, pos_legend, limits = NULL) {
        p_events <- ggplot(df_plot, aes(x = lab, y = val, fill = stat)) +
            # geom_violin(width = 1, alpha = 0.2) +
            geom_split_violin(width = 1.5, alpha = 0.2) +
            geom_boxplot(width = 0.2) +
            coord_flip() +
            labs(fill = "") +
            ylab(lab) +
            scale_fill_manual(
                values = c(
                    "WGD=0" = "magenta4",
                    "WGD=1" = "darkorange"
                ),
                breaks = c("WGD=0", "WGD=1"),
                labels = c("non-WGD cells", "WGD cells")
            ) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 40), legend.position = pos_legend, legend.title = element_text(size = 0), legend.text = element_text(size = 20), axis.title.y = element_text(size = 0), axis.text.y = element_text(size = 0)) +
            guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        if (!is.null(limits)) {
            p_events <- p_events + ylim(limits)
        }
        return(p_events)
    }
    #--------------------Plot the clonal development for all simulations
    filename <- paste0(plotname, ".jpeg")
    jpeg(filename, width = 1000, height = 1100)
    df_stat_sims_all$lab <- factor(df_stat_sims_all$lab, levels = var_labs)
    #   Plot age of MRCA
    df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat == "MRCA_age"), ]
    p_MRCA <- plot_MRCA_age(df_plot, lab = "Age of MRCA")
    #   Plot count of events
    df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("clonal_misseg_loss", "clonal_misseg_gain", "subclonal_misseg_loss", "subclonal_misseg_gain")), ]
    df_plot$stat[which(df_plot$stat == "clonal_misseg_gain")] <- "clonal_gain"
    df_plot$stat[which(df_plot$stat == "clonal_misseg_loss")] <- "clonal_loss"
    df_plot$stat[which(df_plot$stat == "subclonal_misseg_gain")] <- "subclonal_gain"
    df_plot$stat[which(df_plot$stat == "subclonal_misseg_loss")] <- "subclonal_loss"
    p_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts", var_labs = var_labs, pos_legend = c(0.5, 0.11))
    #   Group subplots into a single plot
    p <- grid.arrange(p_MRCA, p_events, widths = c(1.5, 1), nrow = 1)
    #   Print plot
    print(p)
    dev.off()
    #------------------Plot the clonal development for only viable cells
    #-------------------------------------------------in all simulations
    filename <- paste0(plotname, "_VIABLE.jpeg")
    jpeg(filename, width = 1000, height = 1100)
    df_stat_sims_all$lab <- factor(df_stat_sims_all$lab, levels = var_labs)
    #   Plot age of MRCA
    df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat == "MRCA_age"), ]
    p_MRCA <- plot_MRCA_age(df_plot, lab = "Age of MRCA")
    #   Plot count of events
    df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_clonal_misseg_loss", "viable_cells_clonal_misseg_gain", "viable_cells_subclonal_misseg_loss", "viable_cells_subclonal_misseg_gain")), ]
    df_plot$stat[which(df_plot$stat == "viable_cells_clonal_misseg_gain")] <- "clonal_gain"
    df_plot$stat[which(df_plot$stat == "viable_cells_clonal_misseg_loss")] <- "clonal_loss"
    df_plot$stat[which(df_plot$stat == "viable_cells_subclonal_misseg_gain")] <- "subclonal_gain"
    df_plot$stat[which(df_plot$stat == "viable_cells_subclonal_misseg_loss")] <- "subclonal_loss"
    p_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts", var_labs = var_labs, pos_legend = c(0.5, 0.11))
    #   Group subplots into a single plot
    p <- grid.arrange(p_MRCA, p_events, widths = c(1.5, 1), nrow = 1)
    #   Print plot
    print(p)
    dev.off()
    #---------------Plot the clonal development classified by WGD status
    if (plot_WGD == TRUE) {
        filename <- paste0(plotname, "_by_WGD_status_VIABLE.jpeg")
        jpeg(filename, width = 2000, height = 1100)
        clonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_clonal_misseg_gain", "viable_cells_clonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_clonal_misseg_gain", "viable_cells_clonal_misseg_loss"))]))
        subclonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_subclonal_misseg_gain", "viable_cells_subclonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_subclonal_misseg_gain", "viable_cells_subclonal_misseg_loss"))]))
        #   Plot age of MRCA
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat == "MRCA_age"), ]
        p_MRCA <- plot_MRCA_age(df_plot, lab = "Age of MRCA")
        #   Plot fitness
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_WGD=0_fitness", "viable_cells_WGD=1_fitness")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=0_fitness")] <- "WGD=0"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=1_fitness")] <- "WGD=1"
        p_fitness <- plot_fitness_by_WGD_status(df_plot, lab = "Fitness", pos_legend = c(0.25, 0.11))
        #   Plot count of events in non-WGD cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_WGD=0_clonal_misseg_loss", "viable_cells_WGD=0_clonal_misseg_gain", "viable_cells_WGD=0_subclonal_misseg_loss", "viable_cells_WGD=0_subclonal_misseg_gain")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=0_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=0_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=0_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=0_subclonal_misseg_gain")] <- "subclonal_gain"
        p_diploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (non-WGD)", pos_legend = c(0.5, 0.11), limits = clonal_limits)
        #   Plot count of events in tetraploid cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_WGD=1_clonal_misseg_loss", "viable_cells_WGD=1_clonal_misseg_gain", "viable_cells_WGD=1_subclonal_misseg_loss", "viable_cells_WGD=1_subclonal_misseg_gain")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=1_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=1_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=1_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_WGD=1_subclonal_misseg_gain")] <- "subclonal_gain"
        p_tetraploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (WGD)", pos_legend = "none", limits = clonal_limits)
        #   Group subplots into a single plot
        p <- grid.arrange(p_MRCA, p_fitness, p_diploid_events, p_tetraploid_events, widths = c(1.5, 1, 1, 1), nrow = 1)
        #   Print plot
        print(p)
        dev.off()
    }
    #-------------------Plot the clonal development classified by ploidy
    if (plot_WGD == TRUE) {
        #---------------------------------------------Plot for all cells
        filename <- paste0(plotname, "_by_ploidy.jpeg")
        jpeg(filename, width = 2000, height = 1100)
        clonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("clonal_misseg_gain", "clonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("clonal_misseg_gain", "clonal_misseg_loss"))]))
        subclonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("subclonal_misseg_gain", "subclonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("subclonal_misseg_gain", "subclonal_misseg_loss"))]))
        #   Plot age of MRCA
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat == "MRCA_age"), ]
        p_MRCA <- plot_MRCA_age(df_plot, lab = "Age of MRCA")
        #   Plot fitness
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("ploidy=2_fitness", "ploidy=4_fitness")), ]
        df_plot$stat[which(df_plot$stat == "ploidy=2_fitness")] <- "ploidy=2"
        df_plot$stat[which(df_plot$stat == "ploidy=4_fitness")] <- "ploidy=4"
        p_fitness <- plot_fitness_by_ploidy(df_plot, lab = "Fitness", pos_legend = c(0.25, 0.11))
        #   Plot count of events in diploid cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("ploidy=2_clonal_misseg_loss", "ploidy=2_clonal_misseg_gain", "ploidy=2_subclonal_misseg_loss", "ploidy=2_subclonal_misseg_gain")), ]
        df_plot$stat[which(df_plot$stat == "ploidy=2_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "ploidy=2_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "ploidy=2_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "ploidy=2_subclonal_misseg_gain")] <- "subclonal_gain"
        p_diploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (ploidy=2)", pos_legend = c(0.5, 0.11), limits = clonal_limits)
        #   Plot count of events in tetraploid cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("ploidy=4_clonal_misseg_loss", "ploidy=4_clonal_misseg_gain", "ploidy=4_subclonal_misseg_loss", "ploidy=4_subclonal_misseg_gain")), ]
        df_plot$stat[which(df_plot$stat == "ploidy=4_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "ploidy=4_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "ploidy=4_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "ploidy=4_subclonal_misseg_gain")] <- "subclonal_gain"
        p_tetraploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (ploidy=4)", pos_legend = "none", limits = clonal_limits)
        #   Group subplots into a single plot
        p <- grid.arrange(p_MRCA, p_fitness, p_diploid_events, p_tetraploid_events, widths = c(1.5, 1, 1, 1), nrow = 1)
        #   Print plot
        print(p)
        dev.off()
        #-------------------------------------Plot for only viable cells
        filename <- paste0(plotname, "_by_ploidy_VIABLE.jpeg")
        jpeg(filename, width = 2000, height = 1100)
        clonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_clonal_misseg_gain", "viable_cells_clonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_clonal_misseg_gain", "viable_cells_clonal_misseg_loss"))]))
        subclonal_limits <- c(min(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_subclonal_misseg_gain", "viable_cells_subclonal_misseg_loss"))]), max(df_stat_sims_all$val[which(df_stat_sims_all$stat %in% c("viable_cells_subclonal_misseg_gain", "viable_cells_subclonal_misseg_loss"))]))
        #   Plot age of MRCA
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat == "MRCA_age"), ]
        p_MRCA <- plot_MRCA_age(df_plot, lab = "Age of MRCA")
        #   Plot fitness
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_ploidy=2_fitness", "viable_cells_ploidy=4_fitness")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=2_fitness")] <- "ploidy=2"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=4_fitness")] <- "ploidy=4"
        p_fitness <- plot_fitness_by_ploidy(df_plot, lab = "Fitness", pos_legend = c(0.25, 0.11))
        #   Plot count of events in diploid cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_ploidy=2_clonal_misseg_loss", "viable_cells_ploidy=2_clonal_misseg_gain", "viable_cells_ploidy=2_subclonal_misseg_loss", "viable_cells_ploidy=2_subclonal_misseg_gain")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=2_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=2_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=2_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=2_subclonal_misseg_gain")] <- "subclonal_gain"
        p_diploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (ploidy=2)", pos_legend = c(0.5, 0.11), limits = clonal_limits)
        #   Plot count of events in tetraploid cells
        df_plot <- df_stat_sims_all[which(df_stat_sims_all$stat %in% c("viable_cells_ploidy=4_clonal_misseg_loss", "viable_cells_ploidy=4_clonal_misseg_gain", "viable_cells_ploidy=4_subclonal_misseg_loss", "viable_cells_ploidy=4_subclonal_misseg_gain")), ]
        if (length(which(is.nan(df_plot$val))) > 0) df_plot <- df_plot[-which(is.nan(df_plot$val)), ]
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=4_clonal_misseg_loss")] <- "clonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=4_clonal_misseg_gain")] <- "clonal_gain"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=4_subclonal_misseg_loss")] <- "subclonal_loss"
        df_plot$stat[which(df_plot$stat == "viable_cells_ploidy=4_subclonal_misseg_gain")] <- "subclonal_gain"
        p_tetraploid_events <- plot_clonal_and_subclonal_events(df_plot, lab = "Event counts (ploidy=4)", pos_legend = "none", limits = clonal_limits)
        #   Group subplots into a single plot
        p <- grid.arrange(p_MRCA, p_fitness, p_diploid_events, p_tetraploid_events, widths = c(1.5, 1, 1, 1), nrow = 1)
        #   Print plot
        print(p)
        dev.off()
    }
}

simulation_statistics_one_simulation <- function(filename, var1_name, var2_name, var1, var2, sim, plot_WGD) {
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
    evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
    evolution_traj_clonal_ID <- simulation$clonal_evolution$evolution_traj_clonal_ID
    T_current <- simulation$clonal_evolution$T_current

    phylogeny_origin <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_origin
    phylogeny_deathtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_deathtime
    #-------------------Function to find clonal ancestry of given clones
    find_clonal_ancestry <- function(list_subclonal_ancestry) {
        if (length(list_subclonal_ancestry) == 0) {
            clonal_ancestry <- c()
        } else if (length(list_subclonal_ancestry) == 1) {
            clonal_ancestry <- list_subclonal_ancestry[[1]]
        } else {
            clonal_ancestry <- list_subclonal_ancestry[[1]]
            for (i in 2:length(list_subclonal_ancestry)) {
                clonal_ancestry <- intersect(clonal_ancestry, list_subclonal_ancestry[[i]])
            }
        }
        return(clonal_ancestry)
    }
    #-------------------------------------------------Find unique clones
    table_clone <- as.data.frame(table(ID = simulation$sample$sample_clone_ID))
    table_clone$ID_unique <- 0
    table_clone$ID_unique[1] <- 1
    ID_unique <- 1
    if (nrow(table_clone) > 1) {
        for (clone_new in 2:nrow(table_clone)) {
            ID_unique_new <- 0
            for (clone_old in 1:(clone_new - 1)) {
                #   Check if selection rate is same
                if (genotype_list_selection_rate[clone_new] != genotype_list_selection_rate[clone_old]) next
                #   Check if CN profile is same
                if (any(genotype_list_ploidy_chrom[[clone_new]] != genotype_list_ploidy_chrom[[clone_old]])) next
                tmp <- 1
                for (chrom in 1:length(genotype_list_ploidy_chrom[[clone_new]])) {
                    if (genotype_list_ploidy_chrom[[clone_new]][chrom] <= 0) next
                    for (strand in 1:genotype_list_ploidy_chrom[[clone_new]][chrom]) {
                        if (!setequal(genotype_list_ploidy_block[[clone_new]][[chrom]][[strand]], genotype_list_ploidy_block[[clone_old]][[chrom]][[strand]])) tmp <- 0
                    }
                }
                if (tmp == 0) next
                ID_unique_new <- table_clone$ID_unique[clone_old]
            }
            if (ID_unique_new == 0) {
                ID_unique <- ID_unique + 1
                ID_unique_new <- ID_unique
            }
            table_clone$ID_unique[clone_new] <- ID_unique_new
        }
    }
    table_clone_unique <- data.frame(ID_unique = 1:ID_unique)
    table_clone_unique$Freq <- 0
    for (ID_unique in 1:nrow(table_clone_unique)) {
        table_clone_unique$Freq[ID_unique] <- sum(table_clone$Freq[which(table_clone$ID_unique == ID_unique)])
    }
    Clone_ID <- as.numeric(as.vector(table_clone$ID))
    Clone_ID_unique <- as.numeric(as.vector(table_clone_unique$ID_unique))
    #---------------------------------------Find ancestry of every clone
    subclonal_ancestry <- vector("list", length(Clone_ID))
    for (i in 1:length(Clone_ID)) {
        ancestry <- Clone_ID[i]
        while (ancestry[1] != 0) ancestry <- c(evolution_origin[ancestry[1]], ancestry)
        subclonal_ancestry[[i]] <- ancestry
    }
    #-----------------------------------------Find ploidy for each clone
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
    N_chromosomes <- length(genotype_list_ploidy_chrom[[1]])
    vec_CN_block_no <- rep(0, N_chromosomes)
    for (chrom in 1:N_chromosomes) {
        vec_CN_block_no[chrom] <- length(genotype_list_ploidy_block[[1]][[chrom]][[1]])
    }
    #   Find ploidy and selection rate for each unique clone and its mother clone
    table_clone$ploidy <- 0
    table_clone$mother_ploidy <- 0
    for (i in 1:nrow(table_clone)) {
        table_clone$ploidy[i] <- compute_ploidy(vec_CN_block_no, genotype_list_ploidy_chrom[[Clone_ID[i]]], genotype_list_ploidy_block[[Clone_ID[i]]])
        if (evolution_origin[Clone_ID[i]] <= 0) {
            mother_clone <- Clone_ID[i]
        } else {
            mother_clone <- evolution_origin[Clone_ID[i]]
        }
        table_clone$mother_ploidy[i] <- compute_ploidy(vec_CN_block_no, genotype_list_ploidy_chrom[[mother_clone]], genotype_list_ploidy_block[[mother_clone]])
    }
    table_clone$rounded_ploidy <- round(table_clone$ploidy)
    table_clone$rounded_mother_ploidy <- round(table_clone$mother_ploidy)
    ploidy_unique <- unique(table_clone$rounded_ploidy)
    # #-------------------------------------Find WGD status for each clone
    # table_clone$WGD_status <- 0
    # for (i in 1:nrow(table_clone)) {
    #     WGD_status <- 0
    #     clone <- Clone_ID[i]
    #     while (clone > 0) {
    #         if (length(evolution_genotype_changes[[clone]]) > 0) {
    #             for (j in 1:length(evolution_genotype_changes[[clone]])) {
    #                 if (evolution_genotype_changes[[clone]][[j]][1] == "whole-genome-duplication") WGD_status <- 1
    #             }
    #         }
    #         clone <- evolution_origin[clone]
    #     }
    #     table_clone$WGD_status[i] <- WGD_status
    # }
    #-------------------------------------Find WGD status for each clone
    find_WGD_status <- function(subclone_ancestry) {
        WGD_status <- 0
        j <- 0
        while ((WGD_status == 0) & (j < length(subclone_ancestry))) {
            j <- j + 1
            clone_node <- subclone_ancestry[j]
            if (clone_node <= 0) next
            events <- evolution_genotype_changes[[clone_node]]
            if (length(events) > 0) {
                for (event in 1:length(events)) {
                    if (events[[event]][1] == "whole-genome-duplication") {
                        WGD_status <- 1
                    }
                }
            }
        }
        return(WGD_status)
    }
    #   Compute WGD status for each unique clone
    table_clone$WGD_status <- 0
    for (i in 1:nrow(table_clone)) {
        subclone_ancestry <- subclonal_ancestry[[i]]
        table_clone$WGD_status[i] <- find_WGD_status(subclone_ancestry)
    }
    #--------------------------------Statistics: status of WGD clonality
    if (plot_WGD == TRUE) {
        if (max(table_clone$WGD_status == 1) & min(table_clone$WGD_status) == 1) {
            flag_clonal_WGD <- 1
            flag_subclonal_WGD <- 0
        } else if ((max(table_clone$WGD_status == 1) == 1) & (min(table_clone$WGD_status == 1) == 0)) {
            flag_clonal_WGD <- 0
            flag_subclonal_WGD <- 1
        } else {
            flag_clonal_WGD <- 0
            flag_subclonal_WGD <- 0
        }
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "flag_clonal_WGD", flag_clonal_WGD)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "flag_subclonal_WGD", flag_subclonal_WGD)
    }
    #---------------------------------------Statistics: total cell count
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "all_cell_count", sum(table_clone$Freq))
    #-----------------------------Statistics: cell count for each ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_cell_count"), sum(table_clone$Freq[which(table_clone$rounded_ploidy == ploidy)]))
    }
    #-------------------------Statistics: cell count for each WGD status
    if (plot_WGD == TRUE) {
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "WGD=0_cell_count", sum(table_clone$Freq[which(table_clone$WGD_status == 0)]))
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "WGD=1_cell_count", sum(table_clone$Freq[which(table_clone$WGD_status == 1)]))
    }
    #-----------------------Statistics: CCF of 5 biggest clones + others
    vec_CCF_unique <- 100 * table_clone_unique$Freq / sum(table_clone_unique$Freq)
    if (length(vec_CCF_unique) <= 5) {
        vec_CCF <- c(sort(vec_CCF_unique, decreasing = TRUE), rep(0, len = (6 - length(vec_CCF_unique))))
    } else {
        tmp <- sort(vec_CCF_unique, decreasing = TRUE)
        vec_CCF <- c(tmp[1:5], sum(tmp[6:length(tmp)]))
    }
    for (i in 1:5) {
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("CCF_", i), vec_CCF[i])
    }
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "CCF_other", vec_CCF[6])
    #--------------------------------------------Statistics: Age of MRCA
    MRCA_age <- 1 - mean(phylogeny_deathtime[which(phylogeny_origin == 0)]) / T_current
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "MRCA_age", MRCA_age)
    #------------------Statistics: fitness with respect to diploid clone
    #-----------------------------------------------for all alive clones
    table_clone$fitness <- genotype_list_selection_rate[Clone_ID] / genotype_list_selection_rate[1]
    fitness <- sum((table_clone$fitness) * (table_clone$Freq)) / sum(table_clone$Freq)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "fitness", fitness)
    #-----------------------Statistics: count of clonal/subclonal events
    #-----------------------------------------------for all alive clones
    find_event_count <- function(ancestry, event_type, event_subtype = NULL) {
        event_count <- 0
        if (length(ancestry) > 0) {
            for (i in 1:length(ancestry)) {
                if (ancestry[i] == 0) next
                if (length(evolution_genotype_changes[[ancestry[i]]]) == 0) next
                for (j in 1:length(evolution_genotype_changes[[ancestry[i]]])) {
                    if (evolution_genotype_changes[[ancestry[i]]][[j]][1] == event_type) {
                        if (!is.null(event_subtype)) {
                            if (strtoi(evolution_genotype_changes[[ancestry[i]]][[j]][4]) == event_subtype) {
                                event_count <- event_count + 1
                            }
                        } else {
                            event_count <- event_count + 1
                        }
                    }
                }
            }
        }
        return(event_count)
    }
    #   Find count of clonal misseg gains and losses
    clonal_ancestry <- find_clonal_ancestry(subclonal_ancestry)
    clonal_count_missegregation_gain <- find_event_count(clonal_ancestry, "missegregation", 1)
    clonal_count_missegregation_loss <- find_event_count(clonal_ancestry, "missegregation", -1)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clonal_misseg_loss", clonal_count_missegregation_loss)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "clonal_misseg_gain", clonal_count_missegregation_gain)
    #   Find average count of subclonal misseg gains and losses
    table_clone$count_missegregation_gain <- 0
    table_clone$count_missegregation_loss <- 0
    for (k in 1:nrow(table_clone)) {
        table_clone$count_missegregation_gain[k] <- find_event_count(subclonal_ancestry[[k]], "missegregation", 1)
        table_clone$count_missegregation_loss[k] <- find_event_count(subclonal_ancestry[[k]], "missegregation", -1)
    }
    subclonal_count_missegregation_gain <- sum((table_clone$count_missegregation_gain) * (table_clone$Freq)) / sum(table_clone$Freq) - clonal_count_missegregation_gain
    subclonal_count_missegregation_loss <- sum((table_clone$count_missegregation_loss) * (table_clone$Freq)) / sum(table_clone$Freq) - clonal_count_missegregation_loss
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "subclonal_misseg_loss", subclonal_count_missegregation_loss)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "subclonal_misseg_gain", subclonal_count_missegregation_gain)
    #   Find count of clonal misseg gains and losses for only viable cells
    mini_table_clone <- table_clone[which(table_clone$fitness > 0), ]
    mini_subclonal_ancestry <- subclonal_ancestry[which(table_clone$fitness > 0)]
    mini_clonal_ancestry <- find_clonal_ancestry(mini_subclonal_ancestry)
    viable_cells_clonal_count_missegregation_gain <- find_event_count(mini_clonal_ancestry, "missegregation", 1)
    viable_cells_clonal_count_missegregation_loss <- find_event_count(mini_clonal_ancestry, "missegregation", -1)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viable_cells_clonal_misseg_loss", viable_cells_clonal_count_missegregation_loss)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viable_cells_clonal_misseg_gain", viable_cells_clonal_count_missegregation_gain)
    #   Find average count of subclonal misseg gains and losses for only viable cells
    viable_cells_subclonal_count_missegregation_gain <- sum((mini_table_clone$count_missegregation_gain) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - viable_cells_clonal_count_missegregation_gain
    viable_cells_subclonal_count_missegregation_loss <- sum((mini_table_clone$count_missegregation_loss) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - viable_cells_clonal_count_missegregation_loss
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viable_cells_subclonal_misseg_loss", viable_cells_subclonal_count_missegregation_loss)
    df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, "viable_cells_subclonal_misseg_gain", viable_cells_subclonal_count_missegregation_gain)
    #-----------------------Statistics: count of clonal/subclonal events
    #-----------------------------------------for clones based on ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        mini_table_clone <- table_clone[which(table_clone$rounded_ploidy == ploidy), ]
        mini_subclonal_ancestry <- subclonal_ancestry[which(table_clone$rounded_ploidy == ploidy)]
        #   Find count of clonal misseg gains and losses
        clonal_ancestry <- find_clonal_ancestry(mini_subclonal_ancestry)
        clonal_count_missegregation_gain <- find_event_count(clonal_ancestry, "missegregation", 1)
        clonal_count_missegregation_loss <- find_event_count(clonal_ancestry, "missegregation", -1)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_clonal_misseg_loss"), clonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_clonal_misseg_gain"), clonal_count_missegregation_gain)
        #   Find average count of subclonal misseg gains and losses
        subclonal_count_missegregation_gain <- sum((mini_table_clone$count_missegregation_gain) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_gain
        subclonal_count_missegregation_loss <- sum((mini_table_clone$count_missegregation_loss) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_loss
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_subclonal_misseg_loss"), subclonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_subclonal_misseg_gain"), subclonal_count_missegregation_gain)
    }
    #-----------------------Statistics: count of clonal/subclonal events
    #-----------------------------for only viable clones based on ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        mini_table_clone <- table_clone[which(table_clone$rounded_ploidy == ploidy & table_clone$fitness > 0), ]
        mini_subclonal_ancestry <- subclonal_ancestry[which(table_clone$rounded_ploidy == ploidy & table_clone$fitness > 0)]
        #   Find count of clonal misseg gains and losses
        clonal_ancestry <- find_clonal_ancestry(mini_subclonal_ancestry)
        clonal_count_missegregation_gain <- find_event_count(clonal_ancestry, "missegregation", 1)
        clonal_count_missegregation_loss <- find_event_count(clonal_ancestry, "missegregation", -1)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_ploidy=", ploidy, "_clonal_misseg_loss"), clonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_ploidy=", ploidy, "_clonal_misseg_gain"), clonal_count_missegregation_gain)
        #   Find average count of subclonal misseg gains and losses
        subclonal_count_missegregation_gain <- sum((mini_table_clone$count_missegregation_gain) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_gain
        subclonal_count_missegregation_loss <- sum((mini_table_clone$count_missegregation_loss) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_loss
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_ploidy=", ploidy, "_subclonal_misseg_loss"), subclonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_ploidy=", ploidy, "_subclonal_misseg_gain"), subclonal_count_missegregation_gain)
    }
    #-----------------------Statistics: count of clonal/subclonal events
    #-------------------------for only viable clones based on WGD status
    WGD_status_unique <- c(0, 1)
    for (i in 1:length(WGD_status_unique)) {
        WGD_status <- WGD_status_unique[i]
        mini_table_clone <- table_clone[which(table_clone$WGD_status == WGD_status & table_clone$fitness > 0), ]
        mini_subclonal_ancestry <- subclonal_ancestry[which(table_clone$WGD_status == WGD_status & table_clone$fitness > 0)]
        #   Find count of clonal misseg gains and losses
        clonal_ancestry <- find_clonal_ancestry(mini_subclonal_ancestry)
        clonal_count_missegregation_gain <- find_event_count(clonal_ancestry, "missegregation", 1)
        clonal_count_missegregation_loss <- find_event_count(clonal_ancestry, "missegregation", -1)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_WGD=", WGD_status, "_clonal_misseg_loss"), clonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_WGD=", WGD_status, "_clonal_misseg_gain"), clonal_count_missegregation_gain)
        #   Find average count of subclonal misseg gains and losses
        subclonal_count_missegregation_gain <- sum((mini_table_clone$count_missegregation_gain) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_gain
        subclonal_count_missegregation_loss <- sum((mini_table_clone$count_missegregation_loss) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq) - clonal_count_missegregation_loss
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_WGD=", WGD_status, "_subclonal_misseg_loss"), subclonal_count_missegregation_loss)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_WGD=", WGD_status, "_subclonal_misseg_gain"), subclonal_count_missegregation_gain)
    }
    #------------------Statistics: fitness with respect to diploid clone
    #-----------------------------------------for clones based on ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        mini_table_clone <- table_clone[which(table_clone$rounded_ploidy == ploidy), ]
        fitness <- sum((mini_table_clone$fitness) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("ploidy=", ploidy, "_fitness"), fitness)
    }
    #------------------Statistics: fitness with respect to diploid clone
    #-----------------------------for only viable clones based on ploidy
    for (i in 1:length(ploidy_unique)) {
        ploidy <- ploidy_unique[i]
        mini_table_clone <- table_clone[which(table_clone$rounded_ploidy == ploidy & table_clone$fitness > 0), ]
        fitness <- sum((mini_table_clone$fitness) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_ploidy=", ploidy, "_fitness"), fitness)
    }
    #------------------Statistics: fitness with respect to diploid clone
    #-------------------------for only viable clones based on WGD status
    for (i in 1:length(WGD_status_unique)) {
        WGD_status <- WGD_status_unique[i]
        mini_table_clone <- table_clone[which(table_clone$WGD_status == WGD_status & table_clone$fitness > 0), ]
        fitness <- sum((mini_table_clone$fitness) * (mini_table_clone$Freq)) / sum(mini_table_clone$Freq)
        df_stat_sim[nrow(df_stat_sim) + 1, ] <- c(var1, var2, sim, paste0("viable_cells_WGD=", WGD_status, "_fitness"), fitness)
    }
    #--------------------------Return the statistics for this simulation
    return(df_stat_sim)
}
