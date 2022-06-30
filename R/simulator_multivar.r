simulator_multivar <- function(model_prefix = "",
                               model_variables_base = list(),
                               var1_name = "",
                               var1_vals = c(),
                               var2_name = "",
                               var2_vals = c(),
                               n_simulations_per_batch = 0,
                               compute_parallel = FALSE,
                               stage_final = 3,
                               n_clones_min = 1,
                               n_clones_max = Inf) {
    # =====================================SIMULATE MATRIX OF STATISTICS
    mat_simulation_statistics <- matrix(list(), nrow = length(var1_vals), ncol = length(var2_vals))
    pb <- txtProgressBar(
        min = 1, max = length(var1_vals) * length(var2_vals),
        style = 3, width = 50, char = "="
    )
    for (row in 1:length(var1_vals)) {
        for (col in 1:length(var2_vals)) {
            setTxtProgressBar(pb, ((row - 1) * length(var2_vals) + col))
            #   Fix parameters for the entry
            model_variables <- model_variables_base
            model_name <- paste(model_prefix, "-", row, "-", col, sep = "")
            var1 <- var1_vals[row]
            if (var1_name == "rate_WGD") {
                model_variables_batch$general_variables$Value[which(model_variables_batch$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var1
            } else {
                if (var1_name == "rate_missegregation") {
                    model_variables_batch$general_variables$Value[which(model_variables_batch$general_variables$Variable == "prob_CN_missegregation")] <- var1
                }
            }
            var2 <- var2_vals[row]
            if (var2_name == "rate_WGD") {
                model_variables_batch$general_variables$Value[which(model_variables_batch$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- var2
            } else {
                if (var2_name == "rate_missegregation") {
                    model_variables_batch$general_variables$Value[which(model_variables_batch$general_variables$Variable == "prob_CN_missegregation")] <- var2
                }
            }
            #   Save parameter set
            model_name <- paste(model_name_prefix, "-", batch, sep = "")
            SAVE_model_variables(model_name = model_name, model_variables = model_variables_batch)
            #   Run simulations for the entry
            n_simulations <- n_simulations_per_batch
            simulation_statistics <- simulator_full_program(
                model = model_name,
                n_simulations = n_simulations,
                stage_final = stage_final,
                n_clones_min = n_clones_min,
                n_clones_max = n_clones_max,
                save_simulation = FALSE,
                save_newick_tree = FALSE,
                save_cn_profile = FALSE,
                # format_cn_profile = "both",
                internal_nodes_cn_info = FALSE,
                model_readcount = TRUE,
                apply_HMM = TRUE,
                report_progress = FALSE,
                compute_parallel = FALSE,
                pseudo_corrected_readcount = TRUE
            )
        }
    }







    return(mat_simulation_statistics)
}
