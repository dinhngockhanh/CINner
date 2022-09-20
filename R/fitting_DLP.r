#' @export
fitting_DLP <- function(model_name,
                        model_variables,
                        copynumber_DLP,
                        list_parameters,
                        n_cores = NULL) {
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }


    fit_ABC_count <- 1000


    #   Define objective function for ABC fitting
    func_ABC <- function(parameters, model_variables) {
        #   Assign parameters in model variables
        model_variables <- assign_paras_DLP(model_variables, parameters)
        #   Make one simulation
        SIMS <- simulator_full_program(
            model = model_variables,
            model_prefix = "",
            n_simulations = 1,
            stage_final = 2,
            save_simulation = FALSE,
            report_progress = FALSE,
            save_cn_profile = FALSE,
            compute_parallel = FALSE,
            output_variables = c(
                "cn_profiles_wide"
            )
        )
        simulation <- SIMS[[1]]
        copynumber_SIMS <- simulation$sample$cn_profiles_wide
        #   Statistics = ???
        ################################################################
        #####################################################  TEMPORARY
        ##############################  what is the statistics comparing
        ############################  copynumber_SIMS & copynumber_DLP ?
        ################################################################
        stat <- 0
        ################################################################
        return(stat)
    }
    ####################################################################
    #########################################################  TEMPORARY
    #################  for now all parameters are set to be ground truth
    ####################################################################
    #   Simulate table of parameters
    parameters_truth <- as.numeric(c(
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_missegregation"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_chrom_arm_missegregation"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_focal_amplification"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_focal_deletion"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_cnloh_interstitial"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "prob_CN_cnloh_terminal"],
        model_variables$general_variables$Value[model_variables$general_variables$Variable == "rate_driver"],
        model_variables$chromosome_arm_library$s_rate,
        model_variables$driver_library$s_rate
    ))
    sim_param <- matrix(0, nrow = fit_ABC_count, ncol = length(list_parameters))
    for (col in 1:length(list_parameters)) {
        sim_param[, col] <- parameters_truth[col]
        # sim_param[, col] <- runif(fit_ABC_count_arms, min = range_arm_s[1], max = range_arm_s[2])
    }
    ####################################################################
    func_ABC(c(0), model_variables)
}

assign_paras_DLP <- function(model_variables, parameters) {
    return(model_variables)
}
