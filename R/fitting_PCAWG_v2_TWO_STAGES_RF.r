#' @export
fitting_PCAWG <- function(model_name,
                          model_variables,
                          copynumber_PCAWG,
                          n_cores = NULL,
                          n_samples = NULL) {
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    #---------------------Check and correct model variables if necessary
    model_variables <- CHECK_model_variables(model_variables)
    #----------------------Get relevant information from model variables
    # selection_model <- model_variables$selection_model$Value[which(model_variables$selection_model$Variable == "selection_model")]
    cn_info <- model_variables$cn_info
    n_chromosomes <- nrow(model_variables$cn_info)
    prob_CN_whole_genome_duplication <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")])
    prob_CN_missegregation <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")])
    prob_CN_chrom_arm_missegregation <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_chrom_arm_missegregation")])
    prob_CN_focal_amplification <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_amplification")])
    prob_CN_focal_deletion <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_deletion")])
    prob_CN_cnloh_interstitial <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_interstitial")])
    prob_CN_cnloh_terminal <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_terminal")])
    rate_driver <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "rate_driver")])
    #---------Convert model variables for fitting individual chromosomes
    prob_CN_whole_genome_duplication <- prob_CN_whole_genome_duplication / n_chromosomes
    prob_CN_missegregation <- prob_CN_missegregation / n_chromosomes
    prob_CN_chrom_arm_missegregation <- prob_CN_chrom_arm_missegregation / n_chromosomes
    prob_CN_focal_amplification <- prob_CN_focal_amplification / n_chromosomes
    prob_CN_focal_deletion <- prob_CN_focal_deletion / n_chromosomes
    prob_CN_cnloh_interstitial <- prob_CN_cnloh_interstitial / n_chromosomes
    prob_CN_cnloh_terminal <- prob_CN_cnloh_terminal / n_chromosomes
    rate_driver <- rate_driver / n_chromosomes
    #   Model variables for fitting arms
    model_variables_arms <- model_variables
    model_variables_arms$selection_model$Value[which(model_variables_arms$selection_model$Variable == "selection_model")] <- "chrom-arm-and-driver-gene-selection-diploid-base"
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- prob_CN_whole_genome_duplication
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_missegregation")] <- prob_CN_missegregation
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_chrom_arm_missegregation")] <- prob_CN_chrom_arm_missegregation
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_focal_amplification")] <- 0
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_focal_deletion")] <- 0
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_cnloh_interstitial")] <- 0
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "prob_CN_cnloh_terminal")] <- 0
    model_variables_arms$general_variables$Value[which(model_variables_arms$general_variables$Variable == "rate_driver")] <- 0
    #   Model variables for fitting genes
    model_variables_genes <- model_variables
    model_variables_genes$selection_model$Value[which(model_variables_genes$selection_model$Variable == "selection_model")] <- "chrom-arm-and-driver-gene-selection-diploid-base"
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- prob_CN_whole_genome_duplication
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_missegregation")] <- prob_CN_missegregation
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_chrom_arm_missegregation")] <- prob_CN_chrom_arm_missegregation
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_focal_amplification")] <- prob_CN_focal_amplification
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_focal_deletion")] <- prob_CN_focal_deletion
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_cnloh_interstitial")] <- prob_CN_cnloh_interstitial
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "prob_CN_cnloh_terminal")] <- prob_CN_cnloh_terminal
    model_variables_genes$general_variables$Value[which(model_variables_genes$general_variables$Variable == "rate_driver")] <- rate_driver
    #---------------------Find the PCAWG gain/loss map for entire genome
    #   Find genome coordinate from one simulation
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    copynumber_sims <- simulator_full_program(
        model = model_variables_genes,
        model_prefix = "TEST",
        n_simulations = 1,
        stage_final = 2,
        save_simulation = FALSE,
        report_progress = FALSE,
        compute_parallel = FALSE,
        output_variables = c(
            "all_sample_genotype",
            "sample_cell_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile"
        )
    )
    simulation <- copynumber_sims[[1]]
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    CNbins_iteration <- sample_genotype_unique_profile[[1]]
    copynumber_coordinates <- CNbins_iteration[, 1:3]
    copynumber_coordinates$width <- copynumber_coordinates$end - copynumber_coordinates$start + 1
    #   Get gain/loss map from PCAWG for this chromosome at bin level
    PCAWG_delta_genome_bins <- gainloss_PCAWG(
        copynumber_PCAWG,
        copynumber_coordinates,
        ploidy_normalization = TRUE,
        use_rbindlist = TRUE
    )
    #   Get gain/loss map from PCAWG for this chromosome at arm level
    PCAWG_delta_genome_arms <- gainloss_PCAWG(
        copynumber_PCAWG,
        copynumber_coordinates,
        ploidy_normalization = TRUE,
        use_rbindlist = TRUE,
        arm_level = TRUE,
        pos_centromeres = cn_info
    )












    #---------------------------------------------------Fitting with ABC
    list_chromosomes <- model_variables_genes$cn_info$Chromosome
    range_arm_s <- c(0.5, 1.5)
    range_gene_s <- c(1, 1.5)






    # fit_ABC_routine <- "abc-rejection"
    # fit_ABC_count <- 3000
    # fit_ABC_tol <- 0.03

    fit_ABC_routine <- "abc-rf"
    fit_ABC_count <- 1000
    # n_samples <- 200
    # fit_ABC_tol <- 0.03





    #--------------------------------Fit individual chromosomes with ABC
    model_variables_best <- model_variables
    if (is.null(n_samples) == TRUE) {
        PCAWG_N_cases <- length(unique(copynumber_PCAWG$donor_unique_id))
    } else {
        PCAWG_N_cases <- n_samples
    }
    # for (i in 1:length(list_chromosomes)) {
    for (i in 10:10) {
        chromosome_target <- list_chromosomes[i]
        #-----------------Tailor the model variables for this chromosome
        #   Model variables for fitting chromosome arms
        model_variables_chrom_arms <- model_variables_arms
        model_variables_chrom_arms$gc_and_mappability <- model_variables_chrom_arms$gc_and_mappability[which(model_variables_chrom_arms$gc_and_mappability$chr == chromosome_target), ]
        model_variables_chrom_arms$cn_info <- model_variables_chrom_arms$cn_info[which(model_variables_chrom_arms$cn_info$Chromosome == chromosome_target), ]
        model_variables_chrom_arms$driver_library <- model_variables_chrom_arms$driver_library[which(model_variables_chrom_arms$driver_library$Chromosome == chromosome_target), ]
        model_variables_chrom_arms$chromosome_arm_library <- model_variables_chrom_arms$chromosome_arm_library[which(model_variables_chrom_arms$chromosome_arm_library$Chromosome == chromosome_target), ]
        model_variables_chrom_arms$initial_cn <- model_variables_chrom_arms$initial_cn[which(model_variables_chrom_arms$initial_cn$Chromosome == chromosome_target), ]
        #   Model variables for fitting chromosome genes
        model_variables_chrom_genes <- model_variables_genes
        model_variables_chrom_genes$gc_and_mappability <- model_variables_chrom_genes$gc_and_mappability[which(model_variables_chrom_genes$gc_and_mappability$chr == chromosome_target), ]
        model_variables_chrom_genes$cn_info <- model_variables_chrom_genes$cn_info[which(model_variables_chrom_genes$cn_info$Chromosome == chromosome_target), ]
        model_variables_chrom_genes$driver_library <- model_variables_chrom_genes$driver_library[which(model_variables_chrom_genes$driver_library$Chromosome == chromosome_target), ]
        model_variables_chrom_genes$chromosome_arm_library <- model_variables_chrom_genes$chromosome_arm_library[which(model_variables_chrom_genes$chromosome_arm_library$Chromosome == chromosome_target), ]
        model_variables_chrom_genes$initial_cn <- model_variables_chrom_genes$initial_cn[which(model_variables_chrom_genes$initial_cn$Chromosome == chromosome_target), ]
        #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #    !!!!!  LATER: CHANGE model_variables_chrom_genes$initial_others  !!!!
        #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #------------------------Find the PCAWG data for this chromosome
        vec_loc <- which(copynumber_coordinates$chr == chromosome_target)
        #   Restrict PCAWG data on arm level
        PCAWG_delta_arms <- list()
        PCAWG_delta_arms$delta_gain <- PCAWG_delta_genome_arms$delta_gain[vec_loc]
        PCAWG_delta_arms$delta_loss <- PCAWG_delta_genome_arms$delta_loss[vec_loc]
        #   Restrict PCAWG data on bin level
        PCAWG_delta_genes <- list()
        PCAWG_delta_genes$delta_gain <- PCAWG_delta_genome_bins$delta_gain[vec_loc]
        PCAWG_delta_genes$delta_loss <- PCAWG_delta_genome_bins$delta_loss[vec_loc]
        #---------------------------------------Prepare data for fitting
        #---Get list of arms and genes to fit for
        list_arms <- model_variables_chrom_genes$chromosome_arm_library$Arm_ID
        list_genes <- model_variables_chrom_genes$driver_library$Gene_ID
        #---Target statistics = gain/loss map from PCAWG
        #   Get PCAWG gain/loss map on arm level
        target_PCAWG_arms <- c(PCAWG_delta_arms$delta_gain, PCAWG_delta_arms$delta_loss)
        #   Get PCAWG gain/loss map on bin level
        target_PCAWG_genes <- c(PCAWG_delta_genes$delta_gain, PCAWG_delta_genes$delta_loss)
        #---Define objective function for ABC fitting
        func_ABC <- function(parameters, model_variables) {
            parameters <- c(parameters, rep(1, len = (length(list_arms) + length(list_genes) - length(parameters))))
            #   Assign selection rates to arms and genes
            i_para <- 0
            if (length(list_arms) > 0) {
                for (arm in 1:length(list_arms)) {
                    i_para <- i_para + 1
                    model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == list_arms[arm])] <- parameters[i_para]
                }
            }
            if (length(list_genes) > 0) {
                for (gene in 1:length(list_genes)) {
                    i_para <- i_para + 1
                    model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == list_genes[gene])] <- parameters[i_para]
                }
            }
            model_variables <- BUILD_driver_library(
                model_variables = model_variables,
                table_arm_selection_rates = model_variables$chromosome_arm_library,
                table_gene_selection_rates = model_variables$driver_library
            )
            #   Make simulations
            SIMS_chromosome <- simulator_full_program(
                model = model_variables,
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
            SIMS_delta <- gainloss_SIMS(SIMS_chromosome, ploidy_normalization = FALSE)
            #   Statistics = gain/loss map from simulations
            stat <- c(SIMS_delta$delta_gain, SIMS_delta$delta_loss)
            return(stat)
        }
        #--------------Fit the chromosome arms to PCAWG copy number data
        #---Simulate table of parameters
        sim_param <- matrix(0, nrow = fit_ABC_count, ncol = length(list_arms))
        for (col in 1:ncol(sim_param)) {
            sim_param[, col] <- runif(fit_ABC_count, min = range_arm_s[1], max = range_arm_s[2])
        }
        #---Fit by ABC-rejection
        if (fit_ABC_routine == "abc-rejection") {
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
        }
        #---Fit by ABC-random forest
        else if (fit_ABC_routine == "abc-rf") {
            # #   Simulate gain/loss map for each parameter set
            # start_time <- Sys.time()
            # cl <- makePSOCKcluster(n_cores)
            # cat("=================================================================================\n")
            # cat("=================================================================================\n")
            # cat("=================================================================================\n")
            # cat(paste("Fitting arm level of chromosome ", chromosome_target, " with ", numCores, " cores, ABC-random forests with ", fit_ABC_count, " simulations\n", sep = ""))
            # model_variables_chrom_arms <<- model_variables_chrom_arms
            # PCAWG_N_cases <<- PCAWG_N_cases
            # sim_param <<- sim_param
            # func_ABC <<- func_ABC
            # clusterExport(cl, varlist = c(
            #     "model_variables_chrom_arms", "PCAWG_N_cases", "sim_param", "hg19_chrlength",
            #     "func_ABC", "BUILD_driver_library", "simulator_full_program", "gainloss_SIMS", "one_simulation", "createCNmatrix", "normalize_cell_ploidy", "calc_state_mode",
            #     "hc2Newick_MODIFIED", "hc2Newick",
            #     "SIMULATOR_VARIABLES_for_simulation",
            #     "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
            #     "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
            #     "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
            #     "SIMULATOR_FULL_PHASE_2_main", "SIMULATOR_FULL_PHASE_3_main",
            #     "get_cn_profile", "p2_cn_profiles_long", "p2_readcount_model", "rbindlist"
            # ))
            # e <- new.env()
            # e$libs <- .libPaths()
            # clusterExport(cl, "libs", envir = e)
            # clusterEvalQ(cl, .libPaths(libs))
            # library(ape)
            # clusterEvalQ(cl = cl, require(ape))
            # pbo <- pboptions(type = "txt")
            # sim_results_list <- pblapply(cl = cl, X = 1:fit_ABC_count, FUN = function(iteration) {
            #     parameters <- sim_param[iteration, ]
            #     stat <- func_ABC(parameters, model_variables_chrom_arms)
            #     return(stat)
            # })
            # stopCluster(cl)
            # sim_stat <- matrix(0, nrow = fit_ABC_count, ncol = length(target_PCAWG_arms))
            # for (row in 1:fit_ABC_count) {
            #     stat <- sim_results_list[[row]]
            #     sim_stat[row, ] <- stat
            # }
            # end_time <- Sys.time()
            # print(end_time - start_time)
            # #   Prepare ingredients for ABC-random forests
            # vec_para_alias <- paste("para", 1:length(list_arms), sep = "")
            # vec_bin_alias <- paste("gainloss_", 1:length(target_PCAWG_arms), sep = "")
            # obs_rf <- data.frame(matrix(target_PCAWG_arms, nrow = 1))
            # colnames(obs_rf) <- vec_bin_alias
            # #   Perform ABC-random forests
            # ABC_RF <- list()
            # for (para in 1:length(list_arms)) {
            #     #   For each parameter...
            #     #   Create training dataset for this parameter
            #     para_ID <- list_arms[para]
            #     cat(paste("Performing ABC random forests for parameter ", para, "/", ncol(sim_param), " (", para_ID, ")\n", sep = ""))
            #     data_rf <- data.frame(sim_param[, para], sim_stat)
            #     colnames(data_rf) <- c(vec_para_alias[para], vec_bin_alias)
            #     #   Perform regression random forest for this parameter
            #     f <- as.formula(paste(vec_para_alias[para], "~.", sep = ""))
            #     model_rf <- regAbcrf(formula = f, data_rf, paral = TRUE, ncores = n_cores)
            #     #   Get posterior statistics for this parameter
            #     post_rf <- predict(model_rf, obs_rf, data_rf, paral = TRUE, ncores = n_cores)
            #     #   Save results for this parameter
            #     ABC_RF[[paste("data_", para_ID, sep = "")]] <- data_rf
            #     ABC_RF[[paste("model_", para_ID, sep = "")]] <- model_rf
            #     ABC_RF[[paste("obs_", para_ID, sep = "")]] <- obs_rf
            #     ABC_RF[[paste("post_", para_ID, sep = "")]] <- post_rf
            # }
            # #   Save results of ABC-random forests
            # filename <- paste(model_name, "_fitting_chr", chromosome_target, "_arms.rda", sep = "")
            # ABC_RF$target_PCAWG <- target_PCAWG_arms
            # ABC_RF$PCAWG_delta <- PCAWG_delta_arms
            # ABC_RF$sim_param <- sim_param
            # ABC_RF$sim_stat <- sim_stat
            # ABC_RF$list_arms <- list_arms
            # ABC_RF$list_genes <- c()
            # ABC_RF$model_variables <- model_variables_chrom_arms
            # ABC_RF$PCAWG_N_cases <- PCAWG_N_cases
            # ABC_RF$chromosome_arm_library <- model_variables_chrom_arms$chromosome_arm_library
            # ABC_RF$driver_library <- c()
            # save(ABC_RF, file = filename)






            filename <- paste(model_name, "_fitting_chr", chromosome_target, "_arms.rda", sep = "")
            load(filename)
            data_rf <- ABC_RF[["data_10p"]]
            model_rf <- ABC_RF[["model_10p"]]
            obs_rf <- ABC_RF[["obs_10p"]]
            post_rf <- ABC_RF[["post_10p"]]
            print(post_rf)
            get_best_para(data_rf, model_rf, obs_rf, post_rf)





            #---Plot the result of fitting thic chromosome to PCAWG
            filename <- paste(model_name, "_fitting_chr", chromosome_target, "_arms.jpeg", sep = "")
            plot_fitting_PCAWG(filename, ABC_RF, fit_ABC_routine)
        }



        next



        #-------------Fit the chromosome genes to PCAWG copy number data
        #---Find best parameters for chromosome arms
        if (fit_ABC_routine == "abc-rejection") {
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
        } else if (fit_ABC_routine == "abc-rf") {
            param_arms <- rep(1, len = length(list_arms))
            for (para in 1:length(list_arms)) {
                param_arms[para] <- ABC_RF[[paste("post_", list_arms[para], sep = "")]]$expectation
            }
        }
        #---Simulate table of parameters
        sim_param <- matrix(0, nrow = fit_ABC_count, ncol = (length(list_arms) + length(list_genes)))
        for (col in 1:ncol(sim_param)) {
            if (col <= length(list_arms)) {
                #   Assign found selection rates for each arm
                sim_param[, col] <- param_arms[col]
            } else {
                #   Simulate parameters for driver genes
                gene_role <- model_variables_chrom_genes$driver_library$Gene_role[col - length(list_arms)]
                if (gene_role == "ONCOGENE") {
                    sim_param[, col] <- runif(fit_ABC_count, min = range_gene_s[1], max = range_gene_s[2])
                } else if (gene_role == "TSG") {
                    sim_param[, col] <- 1 / runif(fit_ABC_count, min = 1 / range_gene_s[2], max = 1 / range_gene_s[1])
                }
            }
        }
        #---Fit by ABC-rejection
        if (fit_ABC_routine == "abc-rejection") {
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
            # ..........................................................
        }
        #---Fit by ABC-random forest
        else if (fit_ABC_routine == "abc-rf") {
            #   Simulate gain/loss map for each parameter set
            start_time <- Sys.time()
            cl <- makePSOCKcluster(n_cores)
            cat("=================================================================================\n")
            cat("=================================================================================\n")
            cat("=================================================================================\n")
            cat(paste("Fitting gene level of chromosome ", chromosome_target, " with ", numCores, " cores, ABC-random forests with ", fit_ABC_count, " simulations\n", sep = ""))
            model_variables_chrom_genes <<- model_variables_chrom_genes
            PCAWG_N_cases <<- PCAWG_N_cases
            sim_param <<- sim_param
            func_ABC <<- func_ABC
            clusterExport(cl, varlist = c(
                "model_variables_chrom_genes", "PCAWG_N_cases", "sim_param", "hg19_chrlength",
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
            pbo <- pboptions(type = "txt")
            sim_results_list <- pblapply(cl = cl, X = 1:fit_ABC_count, FUN = function(iteration) {
                parameters <- sim_param[iteration, ]
                stat <- func_ABC(parameters, model_variables_chrom_genes)
                return(stat)
            })
            stopCluster(cl)
            sim_stat <- matrix(0, nrow = fit_ABC_count, ncol = length(target_PCAWG_genes))
            for (row in 1:fit_ABC_count) {
                stat <- sim_results_list[[row]]
                sim_stat[row, ] <- stat
            }
            end_time <- Sys.time()
            print(end_time - start_time)
            #   Prepare ingredients for ABC-random forests
            vec_para_alias <- paste("para", 1:length(list_genes), sep = "")
            vec_bin_alias <- paste("gainloss_", 1:length(target_PCAWG_genes), sep = "")
            obs_rf <- data.frame(matrix(target_PCAWG_genes, nrow = 1))
            colnames(obs_rf) <- vec_bin_alias
            #   Perform ABC-random forests
            if (length(list_genes) > 0) {
                for (para in 1:length(list_genes)) {
                    #   For each parameter...
                    #   Create training dataset for this parameter
                    para_ID <- list_genes[para]
                    cat(paste("Performing ABC random forests for parameter ", para, "/", length(list_genes), " (", para_ID, ")\n", sep = ""))
                    data_rf <- data.frame(sim_param[, para + length(list_arms)], sim_stat)
                    colnames(data_rf) <- c(vec_para_alias[para], vec_bin_alias)
                    #   Perform regression random forest for this parameter
                    f <- as.formula(paste(vec_para_alias[para], "~.", sep = ""))
                    model_rf <- regAbcrf(formula = f, data_rf, paral = TRUE, ncores = n_cores)
                    #   Get posterior statistics for this parameter
                    post_rf <- predict(model_rf, obs_rf, data_rf, paral = TRUE, ncores = n_cores)
                    #   Save results for this parameter
                    ABC_RF[[paste("data_", para_ID, sep = "")]] <- data_rf
                    ABC_RF[[paste("model_", para_ID, sep = "")]] <- model_rf
                    ABC_RF[[paste("obs_", para_ID, sep = "")]] <- obs_rf
                    ABC_RF[[paste("post_", para_ID, sep = "")]] <- post_rf
                }
            }
            #   Save results of ABC-random forests
            filename <- paste(model_name, "_fitting_chr", chromosome_target, "_genes_and_arms.rda", sep = "")
            ABC_RF$target_PCAWG <- target_PCAWG_genes
            ABC_RF$PCAWG_delta <- PCAWG_delta_genes
            ABC_RF$sim_param <- sim_param
            ABC_RF$sim_stat <- sim_stat
            ABC_RF$list_arms <- list_arms
            ABC_RF$list_genes <- list_genes
            ABC_RF$model_variables <- model_variables_chrom_genes
            ABC_RF$PCAWG_N_cases <- PCAWG_N_cases
            ABC_RF$chromosome_arm_library <- model_variables_chrom_genes$chromosome_arm_library
            ABC_RF$driver_library <- model_variables_chrom_genes$driver_library
            save(ABC_RF, file = filename)
            #---Plot the result of fitting thic chromosome to PCAWG
            filename <- paste(model_name, "_fitting_chr", chromosome_target, "_genes_and_arms.jpeg", sep = "")
            plot_fitting_PCAWG(filename, ABC_RF, fit_ABC_routine)
        }
        #---Save selection rates predicted by ABC-random forests
        if (length(list_arms) > 0) {
            for (arm in 1:length(list_arms)) {
                para_ID <- list_arms[arm]
                post_rf <- ABC_RF[[paste("post_", para_ID, sep = "")]]
                model_variables_best$chromosome_arm_library$s_rate[which(model_variables_best$chromosome_arm_library$Arm_ID == list_arms[arm])] <- post_rf$expectation
            }
        }
        if (length(list_genes) > 0) {
            for (gene in 1:length(list_genes)) {
                para_ID <- list_genes[gene]
                post_rf <- ABC_RF[[paste("post_", para_ID, sep = "")]]
                model_variables_best$driver_library$s_rate[which(model_variables_best$driver_library$Gene_ID == list_genes[gene])] <- post_rf$expectation
            }
        }
        model_variables_best <- BUILD_driver_library(
            model_variables = model_variables_best,
            table_arm_selection_rates = model_variables_best$chromosome_arm_library,
            table_gene_selection_rates = model_variables_best$driver_library
        )
        ##################   OLD: ABC-REJECTION (TO BE RE-INCORPORATED?)
        # #   Perform ABC-rejection
        # ABC <- abc(target = target_PCAWG, param = sim_param, sumstat = sim_stat, tol = fit_ABC_tol, method = "rejection")
        # #   Save results of ABC-rejection
        # filename <- paste(model_name, "_fitting_chr", chromosome_target, ".rda", sep = "")
        # ABC$target_PCAWG <- target_PCAWG
        # ABC$PCAWG_delta <- PCAWG_delta
        # ABC$sim_param <- sim_param
        # ABC$sim_stat <- sim_stat
        # ABC$list_arms <- list_arms
        # ABC$list_genes <- list_genes
        # ABC$model_variables_chrom_genes <- model_variables_chrom_genes
        # # ABC$PCAWG_chromosome <- PCAWG_chromosome
        # ABC$PCAWG_N_cases <- PCAWG_N_cases
        # ABC$chromosome_arm_library <- model_variables_chrom_genes$chromosome_arm_library
        # ABC$driver_library <- model_variables_chrom_genes$driver_library
        # save(ABC, file = filename)
        # #---Plot the result of fitting thic chromosome to PCAWG
        # filename <- paste(model_name, "_fitting_chr", chromosome_target, ".jpeg", sep = "")
        # plot_fitting_PCAWG(filename, ABC, fit_ABC_routine)
        # #---Save selection rates with best Euclidean score
        # best_param <- sim_param[which(ABC$dist == min(ABC$dist)), ]
        # if (length(list_arms) > 0) {
        #     for (arm in 1:length(list_arms)) {
        #         model_variables_best$chromosome_arm_library$s_rate[which(model_variables_best$chromosome_arm_library$Arm_ID == list_arms[arm])] <- best_param[arm]
        #     }
        # }
        # if (length(list_genes) > 0) {
        #     for (gene in 1:length(list_genes)) {
        #         model_variables_best$driver_library$s_rate[which(model_variables_best$driver_library$Gene_ID == list_genes[gene])] <- best_param[gene + length(list_arms)]
        #     }
        # }
        # model_variables_best <- BUILD_driver_library(
        #     model_variables = model_variables_best,
        #     table_arm_selection_rates = model_variables_best$chromosome_arm_library,
        #     table_gene_selection_rates = model_variables_best$driver_library
        # )
    }
    #--------------------Validate fitted parameter set for entire genome
    cat("Performing validation for entire genome\n")
    #---Save fitted parameter set for entire genome
    filename <- paste(model_name, "_fitting_genome.rda", sep = "")
    save(HGSOC_fitting_genome, file = filename)
    #---Compare simulations against PCAWG for entire genome
    copynumber_sims <- simulator_full_program(
        model = model_variables_best,
        model_prefix = paste(model_name, "_fitting_genome", sep = ""),
        n_simulations = PCAWG_N_cases,
        stage_final = 2,
        save_simulation = TRUE,
        report_progress = TRUE,
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

get_best_para <- function(data_rf,
                          model_rf,
                          obs_rf,
                          post_rf) {
    findweights <- getFromNamespace("findweights", "abcrf")

    object <- model_rf
    obs <- obs_rf
    training <- data_rf
    x <- obs
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[1]
    mf$formula <- object$formula
    mf$data <- training
    training <- mf$data
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    resp <- model.response(mf)
    obj <- object$model.rf
    inbag <- matrix(unlist(obj$inbag.counts, use.names = FALSE), ncol = obj$num.trees, byrow = FALSE)
    ncores <- max(detectCores() - 1, 1)
    obj[["origNodes"]] <- predict(obj, training, predict.all = TRUE, num.threads = max(detectCores() - 1, 1))$predictions
    obj[["origObs"]] <- model.response(mf)
    origObs <- obj$origObs
    origNodes <- obj$origNodes
    nodes <- predict(obj, x, predict.all = TRUE, num.threads = ncores)$predictions
    if (is.null(dim(nodes))) nodes <- matrix(nodes, nrow = 1)
    ntree <- obj$num.trees
    nobs <- object$model.rf$num.samples
    nnew <- nrow(x)

    ntree <- obj$num.trees
    weights <- findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
    weights.std <- weights / ntree

    print(weights.std)
}

#' @export
plot_fitting_PCAWG <- function(filename,
                               ABC,
                               fit_ABC_routine) {
    #----------------------------------------------Input ABC information
    chromosome_arm_library <- ABC$chromosome_arm_library
    driver_library <- ABC$driver_library
    PCAWG_delta <- ABC$PCAWG_delta
    param <- ABC$sim_param
    list_arms <- ABC$list_arms
    list_genes <- ABC$list_genes
    model_variables <- ABC$model_variables
    # PCAWG_chromosome <- ABC$PCAWG_chromosome
    PCAWG_N_cases <- ABC$PCAWG_N_cases
    #   Prepare posterior distributions and best parameter set
    if (fit_ABC_routine == "abc-rejection") {
        #   Prepare posterior distributions
        posterior <- ABC$unadj.values
        #   Prepare best parameter set
        dist <- ABC$dist
        best_param <- param[which(dist == min(dist)), ]
    } else if (fit_ABC_routine == "abc-rf") {
        #   Prepare best parameter set
        best_param <- rep(0, len = (length(list_arms) + length(list_genes)))
        for (para in 1:length(best_param)) {
            if (para <= length(list_arms)) {
                para_ID <- list_arms[para]
            } else {
                para_ID <- list_genes[para - length(list_arms)]
            }
            post_rf <- ABC[[paste("post_", para_ID, sep = "")]]
            best_param[para] <- post_rf$expectation
        }
    }
    #----------------------------Make simulations for best parameter set
    #   Assign selection rates to arms and genes
    i_para <- 0
    if (length(list_arms) > 0) {
        for (arm in 1:length(list_arms)) {
            i_para <- i_para + 1
            model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == list_arms[arm])] <- best_param[i_para]
        }
    }
    if (length(list_genes) > 0) {
        for (gene in 1:length(list_genes)) {
            i_para <- i_para + 1
            model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == list_genes[gene])] <- best_param[i_para]
        }
    }
    model_variables <- BUILD_driver_library(
        model_variables = model_variables,
        table_arm_selection_rates = model_variables$chromosome_arm_library,
        table_gene_selection_rates = model_variables$driver_library
    )
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables,
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
    SIMS_delta <- gainloss_SIMS(SIMS_chromosome,
        ploidy_normalization = FALSE
    )
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
    if (length(which(driver_library$Gene_role == "TSG")) > 0) {
        driver_library_TSG <- driver_library[which(driver_library$Gene_role == "TSG"), ]
        for (gene in 1:nrow(driver_library_TSG)) {
            p_sims <- p_sims + annotate("text",
                x = driver_library_TSG$Bin[gene], y = 0,
                label = driver_library_TSG$Gene_ID[gene],
                size = 10,
                angle = "90",
                hjust = 1
            )
        }
    }
    if (length(which(driver_library$Gene_role == "ONCOGENE")) > 0) {
        driver_library_ONCOGENE <- driver_library[which(driver_library$Gene_role == "ONCOGENE"), ]
        for (gene in 1:nrow(driver_library_ONCOGENE)) {
            p_sims <- p_sims + annotate("text",
                x = driver_library_ONCOGENE$Bin[gene], y = 0,
                label = driver_library_ONCOGENE$Gene_ID[gene],
                size = 10,
                angle = "90",
                hjust = 0
            )
        }
    }
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
            #   Plot distributions of selection rates for this arm
            if (fit_ABC_routine == "abc-rejection") {
                posterior_arm <- posterior[, arm]
                p_arm[[arm]] <- ggplot(data.frame(posterior_arm = posterior_arm), aes(x = posterior_arm)) +
                    geom_density(color = "darkblue", fill = "lightblue") +
                    geom_vline(aes(xintercept = 1), color = "blue", size = 1) +
                    xlab("") +
                    ylab("") +
                    ggtitle(paste("Sigma posterior - ", list_arms[arm], sep = "")) +
                    theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
                    theme(text = element_text(size = 20)) +
                    scale_x_continuous(expand = c(0, 0)) +
                    scale_y_continuous(expand = c(0, 0))
            } else if (fit_ABC_routine == "abc-rf") {
                para_ID <- list_arms[arm]
                data_rf <- ABC[[paste("data_", para_ID, sep = "")]]
                model_rf <- ABC[[paste("model_", para_ID, sep = "")]]
                obs_rf <- ABC[[paste("obs_", para_ID, sep = "")]]
                p_arm[[arm]] <- densityPlot_MODIFIED(
                    model_rf,
                    obs_rf,
                    data_rf,
                    protocol = "arm",
                    chosen_para = best_param[arm],
                    color_prior = "lightblue",
                    color_posterior = "darkblue",
                    color_vline = "blue",
                    main = paste("Sigma posterior - ", list_arms[arm], sep = "")
                )
            }
            #   Plot selected selection rate in best parameter set
            # p_arm[[arm]] <- p_arm[[arm]] +
            #     geom_vline(aes(xintercept = best_param[arm]), color = "blue", size = 1, linetype = "dotted")
        }
    }
    #---Plot posterior distribution for genes
    p_genes <- list()
    if (length(list_genes) > 0) {
        for (gene in 1:length(list_genes)) {
            if (fit_ABC_routine == "abc-rejection") {
                gene_ID <- list_genes[gene]
                gene_role <- driver_library$Gene_role[which(driver_library$Gene_ID == gene_ID)]
                posterior_gene <- posterior[, gene + length(list_arms)]
                if (gene_role == "ONCOGENE") {
                    #   Plot distributions of selection rates for this gene
                    p_genes[[gene]] <- ggplot(data.frame(posterior_gene = posterior_gene), aes(x = posterior_gene)) +
                        geom_density(color = "coral", fill = "lightpink") +
                        geom_vline(aes(xintercept = 1), color = "red", size = 1) +
                        xlab("") +
                        ylab("") +
                        ggtitle(paste("Sigma posterior - ", gene_ID, sep = "")) +
                        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
                        theme(text = element_text(size = 20)) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0))
                    #   Plot selected selection rate in best parameter set
                    p_genes[[gene]] <- p_genes[[gene]] +
                        geom_vline(aes(xintercept = best_param[gene + length(list_arms)]), color = "red", size = 1, linetype = "dotted")
                } else if (gene_role == "TSG") {
                    posterior_gene <- 1 / posterior_gene
                    #   Plot distributions of selection rates for this gene
                    p_genes[[gene]] <- ggplot(data.frame(posterior_gene = posterior_gene), aes(x = posterior_gene)) +
                        geom_density(color = "darkolivegreen4", fill = "darkolivegreen2") +
                        geom_vline(aes(xintercept = 1), color = "forestgreen", size = 1) +
                        xlab("") +
                        ylab("") +
                        ggtitle(paste("Sigma posterior - ", gene_ID, sep = "")) +
                        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
                        theme(text = element_text(size = 20)) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0))
                    #   Plot selected selection rate in best parameter set
                    p_genes[[gene]] <- p_genes[[gene]] +
                        geom_vline(aes(xintercept = (1 / best_param[gene + length(list_arms)])), color = "forestgreen", size = 1, linetype = "dotted")
                }
            } else if (fit_ABC_routine == "abc-rf") {
                gene_ID <- list_genes[gene]
                gene_role <- driver_library$Gene_role[which(driver_library$Gene_ID == gene_ID)]
                para_ID <- list_genes[gene]
                data_rf <- ABC[[paste("data_", para_ID, sep = "")]]
                model_rf <- ABC[[paste("model_", para_ID, sep = "")]]
                obs_rf <- ABC[[paste("obs_", para_ID, sep = "")]]
                if (gene_role == "ONCOGENE") {
                    #   Plot distributions of selection rates for this gene
                    p_genes[[gene]] <- densityPlot_MODIFIED(model_rf,
                        obs_rf,
                        data_rf,
                        protocol = gene_role,
                        chosen_para = best_param[gene + length(list_arms)],
                        color_prior = "lightpink",
                        color_posterior = "coral",
                        color_vline = "red",
                        main = paste("Sigma posterior - ", gene_ID, sep = "")
                    )
                    # #   Plot selected selection rate in best parameter set
                    # p_genes[[gene]] <- p_genes[[gene]] +
                    #     geom_vline(aes(xintercept = best_param[gene + length(list_arms)]), color = "red", size = 1, linetype = "dotted")
                } else if (gene_role == "TSG") {
                    #   Plot distributions of selection rates for this gene
                    p_genes[[gene]] <- densityPlot_MODIFIED(model_rf,
                        obs_rf,
                        data_rf,
                        protocol = gene_role,
                        chosen_para = (1 / best_param[gene + length(list_arms)]),
                        color_prior = "darkolivegreen2",
                        color_posterior = "darkolivegreen4",
                        color_vline = "forestgreen",
                        main = paste("Sigma posterior - ", gene_ID, sep = "")
                    )
                    # #   Plot selected selection rate in best parameter set
                    # p_genes[[gene]] <- p_genes[[gene]] +
                    #     geom_vline(aes(xintercept = (1 / best_param[gene + length(list_arms)])), color = "forestgreen", size = 1, linetype = "dotted")
                }
            }
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
