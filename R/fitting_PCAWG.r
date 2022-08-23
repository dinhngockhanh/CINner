#' @export
fitting_PCAWG <- function(model_name,
                          model_variables,
                          copynumber_PCAWG,
                          n_cores = NULL) {
    #---------------------Check and correct model variables if necessary
    model_variables <- CHECK_model_variables(model_variables)
    model_variables_fitting <- model_variables
    #----------------------Get relevant information from model variables
    selection_model <- model_variables$selection_model$Value[which(model_variables$selection_model$Variable == "selection_model")]
    table_cn <- model_variables$cn_info
    n_chromosomes <- nrow(table_cn)
    prob_CN_whole_genome_duplication <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")])
    prob_CN_missegregation <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")])
    prob_CN_chrom_arm_missegregation <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_chrom_arm_missegregation")])
    prob_CN_focal_amplification <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_amplification")])
    prob_CN_focal_deletion <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_deletion")])
    prob_CN_cnloh_interstitial <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_interstitial")])
    prob_CN_cnloh_terminal <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_terminal")])
    rate_driver <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "rate_driver")])
    #---------Convert model variables for fitting individual chromosomes
    #   Change CN/driver event probabilities to one chromosome
    prob_CN_whole_genome_duplication <- prob_CN_whole_genome_duplication / n_chromosomes
    prob_CN_missegregation <- prob_CN_missegregation / n_chromosomes
    prob_CN_chrom_arm_missegregation <- prob_CN_chrom_arm_missegregation / n_chromosomes
    prob_CN_focal_amplification <- prob_CN_focal_amplification / n_chromosomes
    prob_CN_focal_deletion <- prob_CN_focal_deletion / n_chromosomes
    prob_CN_cnloh_interstitial <- prob_CN_cnloh_interstitial / n_chromosomes
    prob_CN_cnloh_terminal <- prob_CN_cnloh_terminal / n_chromosomes
    rate_driver <- rate_driver / n_chromosomes
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- prob_CN_whole_genome_duplication
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_missegregation")] <- prob_CN_missegregation
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_chrom_arm_missegregation")] <- prob_CN_chrom_arm_missegregation
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_focal_amplification")] <- prob_CN_focal_amplification
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_focal_deletion")] <- prob_CN_focal_deletion
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_cnloh_interstitial")] <- prob_CN_cnloh_interstitial
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "prob_CN_cnloh_terminal")] <- prob_CN_cnloh_terminal
    model_variables_fitting$general_variables$Value[which(model_variables_fitting$general_variables$Variable == "rate_driver")] <- rate_driver
    #---------------------------------------------------Fitting with ABC
    list_chromosomes <- model_variables_fitting$cn_info$Chromosome
    range_arm_s <- c(0.5, 1.5)
    range_gene_s <- c(1, 1.5)





    fit_ABC_count <- 1000
    fit_ABC_tol <- 0.1





    #---Fit individual chromosomes with ABC-rejection
    model_variables_best <- model_variables
    # for (i in 1:length(list_chromosomes)) {
    for (i in 1:length(list_chromosomes)) {
        chromosome_target <- list_chromosomes[i]
        #---Tailor the model variables for this chromosome
        model_variables_chrom <- model_variables_fitting
        model_variables_chrom$gc_and_mappability <- model_variables_chrom$gc_and_mappability[which(model_variables_chrom$gc_and_mappability$chr == chromosome_target), ]
        model_variables_chrom$cn_info <- model_variables_chrom$cn_info[which(model_variables_chrom$cn_info$Chromosome == chromosome_target), ]
        model_variables_chrom$driver_library <- model_variables_chrom$driver_library[which(model_variables_chrom$driver_library$Chromosome == chromosome_target), ]
        model_variables_chrom$chromosome_arm_library <- model_variables_chrom$chromosome_arm_library[which(model_variables_chrom$chromosome_arm_library$Chromosome == chromosome_target), ]
        model_variables_chrom$initial_cn <- model_variables_chrom$initial_cn[which(model_variables_chrom$initial_cn$Chromosome == chromosome_target), ]
        #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #    !!!!!  LATER: CHANGE model_variables_chrom$initial_others  !!!!
        #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #---Find the PCAWG data for this chromosome
        #   Dictate the PCAWG coordinates to specifics of this chromosome
        N_bins <- model_variables_chrom$cn_info$Bin_count
        size_CN_block_DNA <- as.numeric(model_variables_chrom$general_variables$Value[which(model_variables_chrom$general_variables$Variable == "size_CN_block_DNA")])
        copynumber_coordinates <- data.frame(chr = chromosome_target, start = ((0:(N_bins - 1)) * size_CN_block_DNA + 1), end = ((1:N_bins) * size_CN_block_DNA))
        #   Restrict the PCAWG data to this chromosome
        PCAWG_chromosome <- copynumber_PCAWG[which(copynumber_PCAWG$chromosome == chromosome_target), ]
        PCAWG_N_cases <- length(unique(PCAWG_chromosome$donor_unique_id))
        #   Get gain/loss map from PCAWG for this chromosome
        PCAWG_delta <- gainloss_PCAWG(PCAWG_chromosome, copynumber_coordinates)
        #---Fit the chromosome to PCAWG copy number data
        #   Get list of arms and genes to fit for
        list_arms <- model_variables_chrom$chromosome_arm_library$Arm_ID
        list_genes <- model_variables_chrom$driver_library$Gene_ID
        #   Target statistics = gain/loss map from PCAWG
        target_PCAWG <- c(PCAWG_delta$delta_gain, PCAWG_delta$delta_loss)
        #   Define objective function for ABC fitting
        func_ABC <- function(parameters) {
            #   Assign selection rates to arms and genes
            i_para <- 0
            if (length(list_arms) > 0) {
                for (arm in 1:length(list_arms)) {
                    i_para <- i_para + 1
                    model_variables_chrom$chromosome_arm_library$s_rate[which(model_variables_chrom$chromosome_arm_library$Arm_ID == list_arms[arm])] <- parameters[i_para]
                }
            }
            if (length(list_genes) > 0) {
                for (gene in 1:length(list_genes)) {
                    i_para <- i_para + 1
                    model_variables_chrom$driver_library$s_rate[which(model_variables_chrom$driver_library$Gene_ID == list_genes[gene])] <- parameters[i_para]
                }
            }
            model_variables_chrom <- BUILD_driver_library(
                model_variables = model_variables_chrom,
                table_arm_selection_rates = model_variables_chrom$chromosome_arm_library,
                table_gene_selection_rates = model_variables_chrom$driver_library
            )
            #   Make simulations
            SIMS_chromosome <- simulator_full_program(
                model = model_variables_chrom,
                model_prefix = "",
                n_simulations = PCAWG_N_cases,
                stage_final = 2,
                save_simulation = FALSE,
                report_progress = FALSE,
                save_cn_profile = FALSE,
                compute_parallel = FALSE,
                output_variables = c(
                    "all_sample_genotype",
                    "sample_cell_ID",
                    "sample_genotype_unique",
                    "sample_genotype_unique_profile"
                )
            )
            #   Get gain/loss map from the simulations
            SIMS_delta <- gainloss_SIMS(SIMS_chromosome, copynumber_coordinates)
            #   Statistics = gain/loss map from simulations
            stat <- c(SIMS_delta$delta_gain, SIMS_delta$delta_loss)
            return(stat)
        }
        #   Simulate table of parameters
        sim_param <- matrix(0, nrow = fit_ABC_count, ncol = (length(list_arms) + length(list_genes)))
        for (col in 1:ncol(sim_param)) {
            if (col <= 2) {
                sim_param[, col] <- runif(fit_ABC_count, min = range_arm_s[1], max = range_arm_s[2])
            } else {
                sim_param[, col] <- runif(fit_ABC_count, min = range_gene_s[1], max = range_gene_s[2])
            }
        }
        #   Simulate gain/loss map for each parameter set
        start_time <- Sys.time()
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        cat(paste("Fitting chromosome ", chromosome_target, " with ", numCores, " cores, ABC-rejection with ", fit_ABC_count, " simulations, tol=", fit_ABC_tol, "\n", sep = ""))
        model_variables_chrom <<- model_variables_chrom
        PCAWG_N_cases <<- PCAWG_N_cases
        sim_param <<- sim_param
        func_ABC <<- func_ABC
        clusterExport(cl, varlist = c(
            "model_variables_chrom", "PCAWG_N_cases", "sim_param", "hg19_chrlength",
            "func_ABC", "BUILD_driver_library", "simulator_full_program", "gainloss_SIMS", "one_simulation", "createCNmatrix", "normalize_cell_ploidy", "calc_state_mode",
            "hc2Newick_MODIFIED", "hc2Newick",
            "SIMULATOR_VARIABLES_for_simulation",
            "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
            "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
            "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
            "SIMULATOR_FULL_PHASE_2_main", "SIMULATOR_FULL_PHASE_3_main",
            "get_cn_profile", "p2_cn_profiles_long", "p2_readcount_model", "rbindlist"
        ))
        e <- new.env()
        e$libs <- .libPaths()
        clusterExport(cl, "libs", envir = e)
        clusterEvalQ(cl, .libPaths(libs))
        library(ape)
        clusterEvalQ(cl = cl, require(ape))



        sim_results_list <- pblapply(cl = cl, X = 1:fit_ABC_count, FUN = function(iteration) {
            parameters <- sim_param[iteration, ]
            stat <- func_ABC(parameters)
            return(stat)
        })
        sim_stat <- matrix(0, nrow = fit_ABC_count, ncol = length(target_PCAWG))
        for (row in 1:fit_ABC_count) {
            #   Statistics = gain/loss map from simulations
            stat <- sim_results_list[[row]]
            sim_stat[row, ] <- stat
        }
        end_time <- Sys.time()
        print(end_time - start_time)
        #   Perform ABC
        ABC <- abc(target = target_PCAWG, param = sim_param, sumstat = sim_stat, tol = fit_ABC_tol, method = "rejection")
        #   Save results of ABC
        filename <- paste(model_name, "_fitting_chr", chromosome_target, ".rda", sep = "")
        ABC$target_PCAWG <- target_PCAWG
        ABC$PCAWG_delta <- PCAWG_delta
        ABC$sim_param <- sim_param
        ABC$sim_stat <- sim_stat
        ABC$list_arms <- list_arms
        ABC$list_genes <- list_genes
        ABC$model_variables_chrom <- model_variables_chrom
        ABC$PCAWG_chromosome <- PCAWG_chromosome
        save(ABC, file = filename)
        #---Plot the result of fitting thic chromosome to PCAWG
        filename <- paste(model_name, "_fitting_chr", chromosome_target, ".jpeg", sep = "")
        plot_fitting_PCAWG(filename, ABC)
        #---Save selection rates with best Euclidean score
        best_param <- sim_param[which(ABC$dist == min(ABC$dist)), ]
        i_para <- 0
        if (length(list_arms) > 0) {
            for (arm in 1:length(list_arms)) {
                i_para <- i_para + 1
                model_variables_best$chromosome_arm_library$s_rate[which(model_variables_best$chromosome_arm_library$Arm_ID == list_arms[arm])] <- best_param[arm]
            }
        }
        if (length(list_genes) > 0) {
            for (gene in 1:length(list_genes)) {
                i_para <- i_para + 1
                model_variables_best$driver_library$s_rate[which(model_variables_best$driver_library$Gene_ID == list_genes[gene])] <- best_param[gene + 2]
            }
        }
        model_variables_best <- BUILD_driver_library(
            model_variables = model_variables_best,
            table_arm_selection_rates = model_variables_best$chromosome_arm_library,
            table_gene_selection_rates = model_variables_best$driver_library
        )
    }
    #--------------------Validate fitted parameter set for entire genome
    #---Save fitted parameter set for entire genome
    filename <- paste(model_name, "_fitting_genome.rda", sep = "")
    save(model_variables_best, file = filename)
    #---Compare simulations against PCAWG for entire genome
    copynumber_sims <- simulator_full_program(
        model = model_variables_best,
        model_prefix = paste(model_name, "_fitting_genome", sep = ""),
        n_simulations = PCAWG_N_cases,
        stage_final = 2,
        save_simulation = TRUE,
        report_progress = FALSE,
        compute_parallel = TRUE,
        output_variables = c(
            "all_sample_genotype",
            "sample_cell_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile"
        )
    )
    filename <- paste(model_name, "_fitting_genome.jpeg", sep = "")
    plot_gainloss(copynumber_sims, copynumber_PCAWG, filename)
}
#' @export
plot_fitting_PCAWG <- function(filename,
                               ABC) {
    #----------------------------------------------Input ABC information
    PCAWG_delta <- ABC$PCAWG_delta
    posterior <- ABC$unadj.values
    param <- ABC$sim_param
    dist <- ABC$dist
    list_arms <- ABC$list_arms
    list_genes <- ABC$list_genes
    model_variables_chrom <- ABC$model_variables_chrom
    PCAWG_chromosome <- ABC$PCAWG_chromosome
    PCAWG_N_cases <- length(unique(PCAWG_chromosome$donor_unique_id))
    #----------------------------Make simulations for best parameter set
    best_param <- param[which(dist == min(dist)), ]
    #   Assign selection rates to arms and genes
    i_para <- 0
    if (length(list_arms) > 0) {
        for (arm in 1:length(list_arms)) {
            i_para <- i_para + 1
            model_variables_chrom$chromosome_arm_library$s_rate[which(model_variables_chrom$chromosome_arm_library$Arm_ID == list_arms[arm])] <- best_param[i_para]
        }
    }
    if (length(list_genes) > 0) {
        for (gene in 1:length(list_genes)) {
            i_para <- i_para + 1
            model_variables_chrom$driver_library$s_rate[which(model_variables_chrom$driver_library$Gene_ID == list_genes[gene])] <- best_param[i_para]
        }
    }
    model_variables_chrom <- BUILD_driver_library(
        model_variables = model_variables_chrom,
        table_arm_selection_rates = model_variables_chrom$chromosome_arm_library,
        table_gene_selection_rates = model_variables_chrom$driver_library
    )
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables_chrom,
        model_prefix = "",
        n_simulations = PCAWG_N_cases,
        stage_final = 2,
        save_simulation = FALSE,
        report_progress = FALSE,
        compute_parallel = TRUE,
        output_variables = c(
            "all_sample_genotype",
            "sample_cell_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile"
        )
    )
    #   Get gain/loss map from the simulations
    SIMS_delta <- gainloss_SIMS(SIMS_chromosome, copynumber_coordinates)
    #-----------Find Spearman correlations between PCAWG and simulations
    stat_gain <- cor.test(SIMS_delta$delta_gain, PCAWG_delta$delta_gain, method = "spearman", exact = FALSE)
    rho_gain <- stat_gain$estimate[["rho"]]
    pval_gain <- stat_gain$p.value
    stat_loss <- cor.test(SIMS_delta$delta_loss, PCAWG_delta$delta_loss, method = "spearman", exact = FALSE)
    rho_loss <- stat_loss$estimate[["rho"]]
    pval_loss <- stat_loss$p.value
    #--------------------------------------Plot results of fitting PCAWG
    #---Plot gain/loss maps for simulations and PCAWG
    df_plot <- data.frame(
        x = 1:length(SIMS_delta$delta_gain),
        gain_sims = SIMS_delta$delta_gain,
        loss_sims = SIMS_delta$delta_loss,
        gain_data = PCAWG_delta$delta_gain,
        loss_data = PCAWG_delta$delta_loss
    )
    #   Plot gain/loss map for simulations
    p_sims <- ggplot(data = df_plot) +
        geom_bar(aes(x = x, y = gain_sims), stat = "identity", colour = "#E34A33", fill = "#E34A33") +
        geom_bar(aes(x = x, y = loss_sims), stat = "identity", colour = "#3182BD", fill = "#3182BD") +
        geom_hline(yintercept = 0, color = "antiquewhite3") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(
            text = element_text(size = 30),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        scale_x_continuous(limits = c(0, nrow(df_plot)), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ggtitle("Simulations") +
        xlab("") +
        ylab("")
    #   Plot gain/loss map for data
    p_data <- ggplot(data = df_plot) +
        geom_bar(aes(x = x, y = gain_data), stat = "identity", colour = "#E34A33", fill = "#E34A33") +
        geom_bar(aes(x = x, y = loss_data), stat = "identity", colour = "#3182BD", fill = "#3182BD") +
        geom_hline(yintercept = 0, color = "antiquewhite3") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(
            text = element_text(size = 30),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        scale_x_continuous(limits = c(0, nrow(df_plot)), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ggtitle("PCAWG") +
        xlab("") +
        ylab("")
    #   Plot Spearman correlation scores
    p_spearman_gain <- ggplot() +
        annotate("text", x = 0, y = 0, size = 12, colour = "#E34A33", label = paste("rho = ", round(rho_gain, 2), sep = "")) +
        theme_void()
    p_spearman_loss <- ggplot() +
        annotate("text", x = 0, y = 0, size = 12, colour = "#3182BD", label = paste("rho = ", round(rho_loss, 2), sep = "")) +
        theme_void()
    #   Arrange gain/loss maps together
    p_gainloss <- grid.arrange(p_sims, p_data, arrangeGrob(p_spearman_gain, p_spearman_loss, nrow = 1), heights = c(5, 5, 1), ncol = 1)
    #---Plot posterior distributions for chromosome arms
    p_arm <- list()
    if (length(list_arms) > 0) {
        for (arm in 1:length(list_arms)) {
            posterior_arm <- posterior[, arm]
            p_arm[[arm]] <- ggplot(data.frame(posterior_arm = posterior_arm), aes(x = posterior_arm)) +
                # geom_histogram() +
                geom_histogram(color = "coral", fill = "lightpink") +
                geom_vline(aes(xintercept = 1), color = "red", size = 1) +
                # geom_vline(aes(xintercept = mean(vec_sigma_posterior)),
                #     color = "red", linetype = "dashed", size = 1
                # ) +
                xlab("Sigma") +
                ggtitle(paste("Posterior - ", list_arms[arm], sep = "")) +
                theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
                theme(text = element_text(size = 20)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0))
        }
    }
    #---Plot posterior distribution for genes
    p_genes <- list()
    if (length(list_genes) > 0) {
        for (gene in 1:length(list_genes)) {
            posterior_gene <- posterior[, gene + 2]
            p_genes[[gene]] <- ggplot(data.frame(posterior_gene = posterior_gene), aes(x = posterior_gene)) +
                # geom_histogram() +
                geom_histogram(color = "darkblue", fill = "lightblue") +
                geom_vline(aes(xintercept = 1), color = "blue", size = 1) +
                # geom_vline(aes(xintercept = mean(vec_sigma_posterior)),
                #     color = "red", linetype = "dashed", size = 1
                # ) +
                xlab("Sigma") +
                ggtitle(paste("Posterior - ", list_genes[gene], sep = "")) +
                theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
                theme(text = element_text(size = 20)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0))
        }
    }
    #---Arrange the mini-plots and print
    #   Define grid arrange
    N_genes <- length(list_genes)
    N_cols <- 2 + ceiling(N_genes / 2)
    N_rows <- 2
    layout <- matrix(NA, nrow = N_rows, ncol = N_cols)
    gs <- list()
    layout[, 1] <- 1
    gs[[1]] <- p_gainloss
    if (length(list_arms) > 0) {
        for (arm in 1:length(list_arms)) {
            row <- arm
            col <- 2
            id <- arm + 1
            layout[row, col] <- id
            gs[[id]] <- p_arm[[arm]]
        }
    }
    if (N_genes > 0) {
        for (gene in 1:N_genes) {
            row <- (gene - 1) %% 2 + 1
            col <- 3 + floor((gene - 1) / 2)
            id <- gene + 3
            layout[row, col] <- id
            gs[[id]] <- p_genes[[gene]]
        }
    }
    #   Print all mini-plots in grid arrangement
    jpeg(filename, width = 2000, height = 1000)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
}

#' @export
gainloss_PCAWG <- function(copynumber_PCAWG,
                           copynumber_coordinates) {
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------------Get gain/loss consensus from PCAWG
    samplelist_data <- unique(copynumber_PCAWG$donor_unique_id)
    CNbins_list_data <- vector("list", length = length(samplelist_data))
    NA_list_data <- vector("list", length = length(samplelist_data))
    CNbin_length <- copynumber_coordinates$end[1] - copynumber_coordinates$start[1] + 1
    for (iteration in 1:length(samplelist_data)) {
        sample_id <- samplelist_data[iteration]
        #   Get CN data for the sample
        cn_profiles_data <- copynumber_PCAWG[which(copynumber_PCAWG$donor_unique_id == sample_id), ]
        #   Get uniform-bin CN data from the data
        CNbins_iteration <- data.frame(chr = copynumber_coordinates$chr, start = copynumber_coordinates$start, end = copynumber_coordinates$end)
        CNbins_iteration$copy <- 0
        CNbins_iteration$state <- 0
        CNbins_iteration$Min <- 0
        CNbins_iteration$Maj <- 0
        CNbins_iteration$cell_id <- sample_id
        NA_iteration <- rep(1, length = nrow(CNbins_iteration))
        for (row in 1:nrow(cn_profiles_data)) {
            chr <- cn_profiles_data$chromosome[row]
            start <- cn_profiles_data$start[row]
            end <- cn_profiles_data$end[row]
            total_cn <- cn_profiles_data$total_cn[row]
            major_cn <- cn_profiles_data$major_cn[row]
            minor_cn <- cn_profiles_data$minor_cn[row]
            bin_start <- floor(start / CNbin_length) * CNbin_length + 1
            bin_end <- ceiling(start / CNbin_length) * CNbin_length
            vec_loc <- which(CNbins_iteration$chr == chr & CNbins_iteration$start >= bin_start & CNbins_iteration$end <= end)
            CNbins_iteration$copy[vec_loc] <- total_cn
            CNbins_iteration$state[vec_loc] <- total_cn
            CNbins_iteration$Min[vec_loc] <- minor_cn
            CNbins_iteration$Maj[vec_loc] <- major_cn
            NA_iteration[vec_loc] <- 0
        }
        #   Add record for this sample to list
        CNbins_list_data[[iteration]] <- CNbins_iteration
        NA_list_data[[iteration]] <- NA_iteration
    }
    CNbins_data <- rbindlist(CNbins_list_data, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(CNbins_data) <- "data.frame"
    #---Transform CN table into copynumber_sims object
    # copynumber_data <- CNbins_list_data[[1]]
    # copynumber_data <- copynumber_data[, 1:3]
    # copynumber_data$width <- copynumber_data$end - copynumber_data$start + 1
    # for (iteration in 1:length(CNbins_list_data)) {
    #     CNbins_iteration <- CNbins_list_data[[iteration]]
    #     cell_id <- CNbins_iteration$cell_id[1]
    #     copynumber_data[[cell_id]] <- CNbins_iteration$copy
    # }
    copynumber_data <- createCNmatrix(CNbins_data,
        field = plotcol, wholegenome = FALSE,
        fillnaplot = fillna, centromere = FALSE
    )
    #---Normalize ploidy of each sample to 2
    copynumber_data <- normalize_cell_ploidy(copynumber_data)
    #---Replace previously NA in each sample with ploidy 2
    for (sample in 1:length(samplelist_data)) {
        sample_id <- samplelist_data[sample]
        NA_iteration <- NA_list_data[[sample]]
        vec_loc <- which(NA_iteration == 1)
        if (length(vec_loc) > 0) {
            tmp <- copynumber_data[sample_id]
            tmp[vec_loc, ] <- 2
            copynumber_data[sample_id] <- tmp
        }
    }
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_data) - 4
    f1_data <- rowSums(copynumber_data[, 5:ncol(copynumber_data)] > cutoff) / n_samples
    f2_data <- -rowSums(copynumber_data[, 5:ncol(copynumber_data)] < cutoff) / n_samples
    attr(f1_data, "names") <- NULL
    attr(f2_data, "names") <- NULL
    #------------------------------Output gain/loss consensus from PCAWG
    output <- list()
    output$delta_gain <- f1_data
    output$delta_loss <- f2_data
    return(output)
}

#' @export
gainloss_SIMS <- function(copynumber_sims,
                          copynumber_coordinates) {
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------Get gain/loss consensus from simulations
    CNbins_list_sims <- vector("list", length = length(copynumber_sims))
    for (iteration in 1:length(copynumber_sims)) {
        simulation <- copynumber_sims[[iteration]]
        all_sample_genotype <- simulation$sample$all_sample_genotype
        sample_cell_ID <- simulation$sample$sample_cell_ID
        # cn_profiles_long <- simulation$cn_profiles_long
        sample_genotype_unique <- simulation$sample$sample_genotype_unique
        sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
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
        #   Add record for this cell to list
        CNbins_iteration <- sample_genotype_unique_profile[[max_genotype]]
        CNbins_iteration$cell_id <- paste("SIMULATION", iteration, "-Library-1-1", sep = "")
        CNbins_list_sims[[iteration]] <- CNbins_iteration
    }
    # CNbins_sims <- rbindlist(CNbins_list_sims, use.names = FALSE, fill = FALSE, idcol = NULL)
    # class(CNbins_sims) <- "data.frame"
    #---Transform CN table into copynumber_sims object
    copynumber_sims <- CNbins_list_sims[[1]]
    copynumber_sims <- copynumber_sims[, 1:3]
    copynumber_sims$width <- copynumber_sims$end - copynumber_sims$start + 1
    for (iteration in 1:length(CNbins_list_sims)) {
        CNbins_iteration <- CNbins_list_sims[[iteration]]
        cell_id <- CNbins_iteration$cell_id[1]
        copynumber_sims[[cell_id]] <- CNbins_iteration$copy
    }
    # copynumber_sims <- createCNmatrix(CNbins_sims,
    #     field = plotcol, wholegenome = FALSE,
    #     fillnaplot = fillna, centromere = FALSE
    # )
    #---Normalize ploidy of each sample to 2
    copynumber_sims <- normalize_cell_ploidy(copynumber_sims)
    #---Get genome coordinates
    copynumber_coordinates <- copynumber_sims[, 1:4]
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_sims) - 4
    f1_sims <- rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] > cutoff) / n_samples
    f2_sims <- -rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] < cutoff) / n_samples
    attr(f1_sims, "names") <- NULL
    attr(f2_sims, "names") <- NULL
    #------------------------------Output gain/loss consensus from PCAWG
    output <- list()
    output$delta_gain <- f1_sims
    output$delta_loss <- f2_sims
    return(output)
}

#' @export
calc_state_mode <- function(states) {
    state_levels <- unique(states)
    state_mode <- state_levels[
        which.max(tabulate(match(states, state_levels)))
    ]
    if (!is.finite(state_mode)) {
        state_mode <- 2
    }
    return(state_mode)
}

#' @export
normalize_cell_ploidy <- function(copynumber) {
    cell_ids <- colnames(copynumber)
    cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

    for (cell_id in cell_ids) {
        state_mode <- calc_state_mode(copynumber[[cell_id]])
        copynumber[[cell_id]] <- as.integer(ceiling(
            copynumber[[cell_id]] / (state_mode / 2)
        ))
        copynumber[[cell_id]][copynumber[[cell_id]] > 11] <- 11
    }
    return(copynumber)
}
