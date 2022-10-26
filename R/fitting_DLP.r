#' @export
fitting_DLP <- function(model_name,
                        model_variables,
                        copynumber_DLP,
                        list_parameters,
                        n_cores = NULL,
                        ABC_method = "rejection",
                        ABC_simcount = 1000,
                        ABC_tol = 0.1,
                        plot_truth = FALSE) {
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    #-----------------Transform CN table from data into Hungarian format
    copynumber_DLP_mat <- data.matrix(copynumber_DLP, rownames.force = NA)
    copynumber_DLP_mat <- t(copynumber_DLP_mat[, c(4:ncol(copynumber_DLP_mat))])
    #------------------------------------------Get list of parameter IDs
    parameter_IDs <- list_parameters$Variable
    #---------------------------------------------Get true parameter set
    parameters_truth <- c()
    if (plot_truth == TRUE) {
        parameters_truth <- rep(0, length(parameter_IDs))
        for (i in 1:length(parameter_IDs)) {
            parameter_ID <- parameter_IDs[i]
            if (parameter_ID %in% model_variables$general_variables$Variable) {
                parameters_truth[i] <- model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)]
            } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
                parameters_truth[i] <- model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)]
            } else if (parameter_ID %in% model_variables$driver_library$Gene_ID) {
                parameters_truth[i] <- model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == parameter_ID)]
            }
        }
        parameters_truth <- as.numeric(parameters_truth)
    }
    #-----------Define function to assign parameters to proper positions
    assign_paras_DLP <- function(model_variables, parameter_IDs, parameters) {
        for (i in 1:length(parameter_IDs)) {
            parameter_ID <- parameter_IDs[i]
            if (parameter_ID %in% model_variables$general_variables$Variable) {
                model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameters[i]
            } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
                model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)] <- parameters[i]
            } else if (parameter_ID %in% model_variables$driver_library$Gene_ID) {
                model_variables$driver_library$s_rate[which(model_variables$driver_library$Gene_ID == parameter_ID)] <- parameters[i]
            }
        }
        return(model_variables)
    }
    #--------------------------Define objective function for ABC fitting
    func_ABC <- function(parameters, parameter_IDs, copynumber_DLP_mat, model_variables) {
        #   Assign parameters in model variables
        model_variables <- assign_paras_DLP(model_variables, parameter_IDs, parameters)
        #   Compute WT and MUT selection rates for driver genes
        model_variables <- BUILD_driver_library(
            model_variables = model_variables,
            table_arm_selection_rates = model_variables$chromosome_arm_library,
            table_gene_selection_rates = model_variables$driver_library
        )
        #   Make one simulation
        SIMS <- simulator_full_program(
            model = model_variables,
            model_prefix = "",
            n_simulations = 1,
            stage_final = 2,
            save_simulation = FALSE, report_progress = FALSE,
            output_variables = c("cn_profiles_wide")
        )
        simulation <- SIMS[[1]]
        copynumber_SIM <- simulation$sample$cn_profiles_wide
        #   Transform CN table from simulation into Hungarian format
        copynumber_SIM_mat <- data.matrix(copynumber_SIM, rownames.force = NA)
        copynumber_SIM_mat <- t(copynumber_SIM_mat[, c(4:ncol(copynumber_SIM_mat))])
        #   Compute Euclidean distance between every pair of cells in data and simulation
        euclidean <- function(a, b) sqrt(rowSums((a - b)^2))
        distance_mat <- matrix(nrow = nrow(copynumber_DLP_mat), ncol = nrow(copynumber_SIM_mat))
        for (row in 1:nrow(copynumber_DLP_mat)) {
            distance_mat[row, ] <- euclidean(rep(copynumber_DLP_mat[row, ], each = nrow(copynumber_SIM_mat)), copynumber_SIM_mat)
        }
        #   Statistics = CN table in simulation, sorted to match cells in data via the Hungarian algorithm
        pairing <- HungarianSolver(distance_mat)$pairs
        DLP_indices <- pairing[, 1]
        SIM_indices <- pairing[, 2]
        if (!all(sort(DLP_indices, index.return = TRUE)$x == sort(DLP_indices, index.return = TRUE)$ix)) {
            SIM_indices <- SIM_indices[sort(DLP_indices, index.return = TRUE)$ix]
            DLP_indices <- sort(DLP_indices, index.return = TRUE)$x
        }
        copynumber_SIM_mat <- copynumber_SIM_mat[SIM_indices, ]
        return(copynumber_SIM_mat)
    }
    #---------------------------------------Simulate table of parameters
    sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    for (col in 1:ncol(sim_param)) {
        sim_param[, col] <- runif(ABC_simcount, min = as.numeric(list_parameters[col, 2]), max = as.numeric(list_parameters[col, 3]))
    }
    #-----------------------------------------------Make reference table
    start_time <- Sys.time()
    cl <- makePSOCKcluster(n_cores)
    cat("Creating simulated CN tables for ABC...\n")
    sim_param <<- sim_param
    parameter_IDs <<- parameter_IDs
    copynumber_DLP_mat <<- copynumber_DLP_mat
    model_variables <<- model_variables
    func_ABC <<- func_ABC
    assign_paras_DLP <<- assign_paras_DLP
    clusterExport(cl, varlist = c(
        "sim_param", "parameter_IDs", "copynumber_DLP_mat", "model_variables", "func_ABC", "assign_paras_DLP",
        "BUILD_driver_library", "simulator_full_program", "one_simulation",
        "SIMULATOR_VARIABLES_for_simulation",
        "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
        "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
        "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
        "SIMULATOR_FULL_PHASE_2_main",
        "SIMULATOR_FULL_PHASE_3_main",
        "get_cn_profile", "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model",
        "HungarianSolver"
    ))
    e <- new.env()
    e$libs <- .libPaths()
    clusterExport(cl, "libs", envir = e)
    clusterEvalQ(cl, .libPaths(libs))
    library(RcppHungarian)
    clusterEvalQ(cl = cl, require(RcppHungarian))
    pbo <- pboptions(type = "txt")
    sim_results_list <- pblapply(cl = cl, X = 1:ABC_simcount, FUN = function(iteration) {
        parameters <- sim_param[iteration, ]
        stat <- func_ABC(parameters, parameter_IDs, copynumber_DLP_mat, model_variables)
        return(stat)
    })
    stopCluster(cl)
    end_time <- Sys.time()
    print(end_time - start_time)
    #   Save the parameters and their statistics
    ABC_input <- list()
    ABC_input$model_name <- model_name
    ABC_input$copynumber_DLP <- copynumber_DLP
    ABC_input$list_parameters <- list_parameters
    ABC_input$sim_param <- sim_param
    ABC_input$sim_results_list <- sim_results_list
    filename <- paste(model_name, "_ABC_input.rda", sep = "")
    save(ABC_input, file = filename)
    #--------------------------------------------------------Perform ABC
    if (ABC_method == "rejection") {
        #---ABC-REJECTION
        #   Simulated statistics = Euclidean distance between sorted CN tables
        sim_stat <- matrix(0, nrow = ABC_simcount, ncol = 1)
        for (row in 1:ABC_simcount) {
            copynumber_SIM_mat <- sim_results_list[[row]]
            sim_stat[row, ] <- sum(sqrt(rowSums((copynumber_SIM_mat - copynumber_DLP_mat)^2)))
        }
        #   Perform ABC-rejection
        ABC_output <- abc(target = 0, param = sim_param, sumstat = sim_stat, tol = ABC_tol, method = "rejection")
        #   Choose best parameter set from posteriors
    }
    #   Save ABC output package
    ABC_output$ABC_method <- ABC_method
    ABC_output$ABC_simcount <- ABC_simcount
    ABC_output$ABC_tol <- ABC_tol
    filename <- paste(model_name, "_ABC_", ABC_method, "_output.rda", sep = "")
    save(ABC_output, file = filename)
    #---------------------------------------------------Plot ABC results
    plot_fitting_DLP(ABC_input, ABC_output, plot_truth, parameters_truth)
}

plot_fitting_DLP <- function(ABC_input,
                             ABC_output,
                             plot_truth = FALSE,
                             parameters_truth = c()) {
    ABC_method <- ABC_output$ABC_method
    IDs <- ABC_input$list_parameters$Variable
    #   Prepare prior and posterior distributions
    if (ABC_method == "rejection") {
        priors <- ABC_input$sim_param
        posteriors <- ABC_output$unadj.values
    }
    #   Plot prior and posterior for each parameter
    for (i in 1:length(IDs)) {
        ID <- IDs[i]
        filename <- paste("ABC_", ABC_method, "_", ID, ".jpeg", sep = "")
        prior <- priors[, i]
        posterior <- posteriors[, i]
        jpeg(filename, width = 2000, height = 1000)
        df_plot <- data.frame(x = c(prior, posterior), dist = c(rep("prior", length(prior)), rep("posterior", length(posterior))))
        p <- ggplot(df_plot, aes(x = x, fill = dist)) +
            geom_density(alpha = .25) +
            scale_fill_manual(values = c("lightblue", "darkblue")) +
            xlab("") +
            ylab("") +
            ggtitle(ID) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
        if (plot_truth == TRUE) {
            p <- p + geom_vline(aes(xintercept = parameters_truth[i]), color = "darkblue", size = 1, linetype = "dashed")
        }
        print(p)
        dev.off()
    }
}
