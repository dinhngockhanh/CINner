#' Produce simulations and outputs files as requested by user
#'
#' @param model_name ...
#' @param cell_lifespan ...
#' @param T_0 ...
#' @param T_end ...
#' @param T_tau_step ...
#' @param Table_sample ...
#' @param table_population_dynamics ...
#' @param Population_end ...
#' @param Max_events ...
#' @param CN_bin_length ...
#' @param prob_CN_whole_genome_duplication ...
#' @param prob_CN_missegregation ...
#' @param prob_CN_chrom_arm_missegregation ...
#' @param prob_CN_focal_amplification ...
#' @param prob_CN_focal_deletion ...
#' @param prob_CN_cnloh_interstitial ...
#' @param prob_CN_cnloh_terminal ...
#' @param prob_CN_focal_amplification_length ...
#' @param prob_CN_focal_deletion_length ...
#' @param prob_CN_cnloh_interstitial_length ...
#' @param prob_CN_cnloh_terminal_length ...
#' @param rate_driver ...
#' @param rate_passenger ...
#' @param selection_model ...
#' @param bound_driver ...
#' @param bound_average_ploidy ...
#' @param bound_homozygosity ...
#' @param SFS_totalsteps ...
#' @param prob_coverage ...
#' @param alpha_coverage ...
#' @param lower_limit_cell_counts ...
#' @param lower_limit_alt_counts ...
#' @param lower_limit_tot_counts ...
#' @param gc ...
#' @param gc_slope ...
#' @param gc_int ...
#' @param sigma1 ...
#' @param num_reads ...
#' @export
BUILD_general_variables <- function(model_name = "MODEL",
                                    cell_lifespan = 4,
                                    T_0 = list(0, "year"),
                                    T_end = list(100, "year"),
                                    T_tau_step = Inf,
                                    Table_sample = data.frame(),
                                    table_population_dynamics = matrix(0, ncol = 2, nrow = 2),
                                    Population_end = Inf,
                                    Max_events = Inf,
                                    CN_bin_length = 500000,
                                    # model_CN_whole_genome_duplication = "per_division",
                                    # model_CN_missegregation = "per_division",
                                    # model_CN_chrom_arm_missegregation = "per_division",
                                    # model_CN_focal_amplification = "per_division",
                                    # model_CN_focal_deletion = "per_division",
                                    # model_CN_cnloh_interstitial = "per_division",
                                    # model_CN_cnloh_terminal = "per_division",
                                    prob_CN_whole_genome_duplication = 0,
                                    prob_CN_missegregation = 0,
                                    prob_CN_chrom_arm_missegregation = 0,
                                    prob_CN_focal_amplification = 0,
                                    prob_CN_focal_deletion = 0,
                                    prob_CN_cnloh_interstitial = 0,
                                    prob_CN_cnloh_terminal = 0,
                                    formula_CN_whole_genome_duplication = "prob_CN_whole_genome_duplication",
                                    formula_CN_missegregation = "prob_CN_missegregation",
                                    formula_CN_chrom_arm_missegregation = "prob_CN_chrom_arm_missegregation",
                                    formula_CN_focal_amplification = "prob_CN_focal_amplification",
                                    formula_CN_focal_deletion = "prob_CN_focal_deletion",
                                    formula_CN_cnloh_interstitial = "prob_CN_cnloh_interstitial",
                                    formula_CN_cnloh_terminal = "prob_CN_cnloh_terminal",
                                    prob_neutral_CN_whole_genome_duplication = 0,
                                    prob_neutral_CN_missegregation = 0,
                                    prob_neutral_CN_chrom_arm_missegregation = 0,
                                    prob_neutral_CN_focal_amplification = 0,
                                    prob_neutral_CN_focal_deletion = 0,
                                    prob_neutral_CN_cnloh_interstitial = 0,
                                    prob_neutral_CN_cnloh_terminal = 0,
                                    model_CN_focal_amplification_length = "geom",
                                    prob_CN_focal_amplification_length = 0.1,
                                    prob_CN_focal_amplification_length_shape_1 = 0,
                                    prob_CN_focal_amplification_length_shape_2 = 0,
                                    model_CN_focal_deletion_length = "geom",
                                    prob_CN_focal_deletion_length_shape_1 = 0,
                                    prob_CN_focal_deletion_length_shape_2 = 0,
                                    prob_CN_focal_deletion_length = 0.1,
                                    model_CN_cnloh_interstitial_length = "geom",
                                    prob_CN_cnloh_interstitial_length = 0.1,
                                    prob_CN_cnloh_interstitial_length_shape_1 = 0,
                                    prob_CN_cnloh_interstitial_length_shape_2 = 0,
                                    model_CN_cnloh_terminal_length = "geom",
                                    prob_CN_cnloh_terminal_length = 0.1,
                                    prob_CN_cnloh_terminal_length_shape_1 = 0,
                                    prob_CN_cnloh_terminal_length_shape_2 = 0,
                                    rate_driver = 0,
                                    rate_passenger = 0,
                                    selection_model = "",
                                    bound_driver = 3,
                                    bound_average_ploidy = 5,
                                    bound_homozygosity = 0,
                                    bound_maximum_CN = Inf,
                                    bound_maximum_CN_normalized = Inf,
                                    SFS_totalsteps = 25,
                                    prob_coverage = 0.05,
                                    alpha_coverage = 0.7,
                                    lower_limit_cell_counts = 0,
                                    lower_limit_alt_counts = 3,
                                    lower_limit_tot_counts = 0,
                                    gc = data.frame(matrix(ncol = 5, nrow = 0)),
                                    gc_slope = 1.2,
                                    gc_int = 0,
                                    sigma1 = 0.1,
                                    num_reads = 1e6) {
    #-----------------------Build model input file for general variables
    columns <- c("Variable", "Value", "Unit", "Note")
    TABLE_VARIABLES <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(TABLE_VARIABLES) <- columns
    N_row <- 0
    #   Set up the start time of simulations
    age_birth <- T_0[[1]]
    age_birth_unit <- T_0[[2]]
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("age_birth", age_birth, age_birth_unit, "Age when simulation starts")
    if (age_birth_unit == "day") {
        T_start_time <- age_birth
    } else {
        if (age_birth_unit == "week") {
            T_start_time <- 7 * age_birth
        } else {
            if (age_birth_unit == "month") {
                T_start_time <- 30 * age_birth
            } else {
                if (age_birth_unit == "year") {
                    T_start_time <- 365 * age_birth
                }
            }
        }
    }
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("T_start_time", T_start_time, "day", "Age when simulation starts (for internal use)")
    #   Set up the end time of simulations
    age_end <- T_end[[1]]
    age_end_unit <- T_end[[2]]
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("age_end", age_end, age_end_unit, "Age when simulation stops")
    if (age_birth_unit != age_end_unit) {
        print("START AND END TIMES DO NOT USE THE SAME UNIT")
    }
    if (age_end_unit == "day") {
        T_end_time <- age_end
    } else {
        if (age_end_unit == "week") {
            T_end_time <- 7 * age_end
        } else {
            if (age_end_unit == "month") {
                T_end_time <- 30 * age_end
            } else {
                if (age_end_unit == "year") {
                    T_end_time <- 365 * age_end
                }
            }
        }
    }
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("T_end_time", T_end_time, "day", "Age when simulation stops (for internal use)")
    #   Set up the time step for tau-leaping algorithm
    N_row <- N_row + 1
    if (T_tau_step == Inf) {
        T_tau_step <- cell_lifespan / 2
    }
    TABLE_VARIABLES[N_row, ] <- c("T_tau_step", T_tau_step, "day", "Time step for tau-leaping algorithm for simulation")
    #   Set up condition to end simulation prematurely based on population size or event count
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("Population_end", Population_end, "cell count", "Condition for ending simulation (Inf if no condition)")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("Max_events", Max_events, "event count", "Condition for ending simulation (Inf if no condition)")
    #   Set up cell lifespan
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("cell_lifespan", cell_lifespan, "day", "Lifespan of one cell")
    #   Set up CN bin width
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("size_CN_block_DNA", CN_bin_length, "bp", "CN bin width")
    #   Set up CN event probabilities
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_whole_genome_duplication", formula_CN_whole_genome_duplication, "", "Formula of probability for a cell division to harbor a WGD event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_whole_genome_duplication", prob_CN_whole_genome_duplication, "", "Probability for a cell division to harbor a WGD event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_missegregation", formula_CN_missegregation, "", "Formula of probability for a cell division to harbor a chromosome mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_missegregation", prob_CN_missegregation, "", "Probability for a cell division to harbor a chromosome mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_chrom_arm_missegregation", formula_CN_chrom_arm_missegregation, "", "Formula of probability for a cell division to harbor a chromosome-arm mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_chrom_arm_missegregation", prob_CN_chrom_arm_missegregation, "", "Probability for a cell division to harbor a chromosome-arm mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_focal_amplification", formula_CN_focal_amplification, "", "Formula of probability for a cell division to harbor a focal amplification event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_amplification", prob_CN_focal_amplification, "", "Probability for a cell division to harbor a focal amplification event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_focal_deletion", formula_CN_focal_deletion, "", "Formula of probability for a cell division to harbor a focal deletion event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_deletion", prob_CN_focal_deletion, "", "Probability for a cell division to harbor a focal deletion event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_cnloh_interstitial", formula_CN_cnloh_interstitial, "", "Formula of probability for a cell division to harbor an interstitial CN-LOH event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_interstitial", prob_CN_cnloh_interstitial, "", "Probability for a cell division to harbor an interstitial CN-LOH event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("formula_CN_cnloh_terminal", formula_CN_cnloh_terminal, "", "Formula of probability for a cell division to harbor a terminal CN-LOH event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_terminal", prob_CN_cnloh_terminal, "", "Probability for a cell division to harbor a terminal CN-LOH event")
    #   Set up neutral CN event probabilities
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_whole_genome_duplication", prob_neutral_CN_whole_genome_duplication, "per cell division", "Probability for a cell division to harbor a neutral WGD event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_missegregation", prob_neutral_CN_missegregation, "per cell division", "Probability for a cell division to harbor a neutral chromosome mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_chrom_arm_missegregation", prob_neutral_CN_chrom_arm_missegregation, "per cell division", "Probability for a cell division to harbor a neutral chromosome-arm mis-segregation event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_focal_amplification", prob_neutral_CN_focal_amplification, "per cell division", "Probability for a cell division to harbor a neutral focal amplification event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_focal_deletion", prob_neutral_CN_focal_deletion, "per cell division", "Probability for a cell division to harbor a neutral focal deletion event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_cnloh_interstitial", prob_neutral_CN_cnloh_interstitial, "per cell division", "Probability for a cell division to harbor a neutral interstitial CN-LOH event")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_neutral_CN_cnloh_terminal", prob_neutral_CN_cnloh_terminal, "per cell division", "Probability for a cell division to harbor a neutral terminal CN-LOH event")
    #   Set up parameters for distributions of lengths of local CN events
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("model_CN_focal_amplification_length", model_CN_focal_amplification_length, "", "Choice of distribution for the block length of a focal amplification event")
    if (model_CN_focal_amplification_length == "geom") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_amplification_length", prob_CN_focal_amplification_length, "", "Geometric parameter for the block length of a focal amplification event")
    } else if (model_CN_focal_amplification_length == "beta") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_amplification_length_shape_1", prob_CN_focal_amplification_length_shape_1, "", "Beta parameter shape 1 for the block length of a focal amplification event")
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_amplification_length_shape_2", prob_CN_focal_amplification_length_shape_2, "", "Beta parameter shape 2 for the block length of a focal amplification event")
    }
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("model_CN_focal_deletion_length", model_CN_focal_deletion_length, "", "Choice of distribution for the block length of a focal deletion event")
    if (model_CN_focal_deletion_length == "geom") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_deletion_length", prob_CN_focal_deletion_length, "", "Geometric parameter for the block length of a focal deletion event")
    } else if (model_CN_focal_deletion_length == "beta") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_deletion_length_shape_1", prob_CN_focal_deletion_length_shape_1, "", "Beta parameter shape 1 for the block length of a focal deletion event")
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_focal_deletion_length_shape_2", prob_CN_focal_deletion_length_shape_2, "", "Beta parameter shape 2 for the block length of a focal deletion event")
    }
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("model_CN_cnloh_interstitial_length", model_CN_cnloh_interstitial_length, "", "Choice of distribution for the block length of an interstitial CN-LOH event")
    if (model_CN_cnloh_interstitial_length == "geom") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_interstitial_length", prob_CN_cnloh_interstitial_length, "", "Geometric parameter for the block length of an interstitial CN-LOH event")
    } else if (model_CN_cnloh_interstitial_length == "beta") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_interstitial_length_shape_1", prob_CN_cnloh_interstitial_length_shape_1, "", "Beta parameter shape 1 for the block length of an interstitial CN-LOH event")
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_interstitial_length_shape_2", prob_CN_cnloh_interstitial_length_shape_2, "", "Beta parameter shape 2 for the block length of an interstitial CN-LOH event")
    }
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("model_CN_cnloh_terminal_length", model_CN_cnloh_terminal_length, "", "Choice of distribution for the block length of a terminal CN-LOH event")
    if (model_CN_cnloh_terminal_length == "geom") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_terminal_length", prob_CN_cnloh_terminal_length, "", "Geometric parameter for the block length of a terminal CN-LOH event")
    } else if (model_CN_cnloh_terminal_length == "beta") {
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_terminal_length_shape_1", prob_CN_cnloh_terminal_length_shape_1, "", "Beta parameter shape 1 for the block length of a terminal CN-LOH event")
        N_row <- N_row + 1
        TABLE_VARIABLES[N_row, ] <- c("prob_CN_cnloh_terminal_length_shape_2", prob_CN_cnloh_terminal_length_shape_2, "", "Beta parameter shape 2 for the block length of a terminal CN-LOH event")
    }
    #   Set up mutation rates
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("rate_driver", rate_driver, "per bp per cell division", "Poisson rate of getting new driver mutations")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("rate_passenger", rate_passenger, "per bp per cell division", "Poisson rate of getting new passenger mutations")
    #   Set up variables for read count model
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("gc_slope", gc_slope, "", "Slope for linear GC model")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("gc_int", gc_int, "", "Intercept for linear GC model")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("sigma1", sigma1, "", "Gamma scale for read depth noise")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("num_reads", num_reads, "", "Number of reads per cell")
    #   Set up variables for sequencing
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("SFS_totalsteps", SFS_totalsteps, "", "Bin count in SFS data")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("prob_coverage", prob_coverage, "", "Mean coverage depth")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("alpha_coverage", alpha_coverage, "", "Alpha parameter for coverage depth")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("lower_limit_cell_counts", lower_limit_cell_counts, "", "Lower limit of cell counts for mutations to be detected")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("lower_limit_alt_counts", lower_limit_alt_counts, "", "Lower limit of alternate read counts for mutations to be detected")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("lower_limit_tot_counts", lower_limit_tot_counts, "", "Lower limit of total read counts for mutations to be detected")
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    #-------------------------Build model input file for selection model
    columns <- c("Variable", "Value", "Unit", "Note")
    TABLE_SELECTION <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(TABLE_SELECTION) <- columns
    N_row <- 0
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("selection_model", selection_model, "", "Choice of selection model")
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("bound_driver", bound_driver, "driver count", "Maximum driver count in viable cells (cells exceeding this will die)")
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("bound_average_ploidy", bound_average_ploidy, "", "Maximum average ploidy across genome (cells exceeding this will die)")
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("bound_maximum_CN", bound_maximum_CN, "", "Maximum CN at any bin (cells exceeding this will die)")
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("bound_maximum_CN_normalized", bound_maximum_CN_normalized, "", "Maximum CN normalized by average ploidy at any bin (cells exceeding this will die)")
    N_row <- N_row + 1
    TABLE_SELECTION[N_row, ] <- c("bound_homozygosity", bound_homozygosity, "", "Maximum number of bins under homozygosity (cells exceeding this will die)")
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    ##########################
    #-------------------Build model input file for chromosome bin counts
    #-------------------------------------------and centromere locations
    vec_chromosome_name <- c(1:22, "X", "Y")
    vec_chromosome_bp <- c(
        248956422, 242193529, 198295559, 190214555,
        181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189, 90338345,
        83257441, 80373285, 58617616, 64444167,
        46709983, 50818468, 156040895, 57227415
    )
    vec_centromere_bp <- c(
        125, 93.3, 91, 50.4,
        48.4, 61, 59.9, 45.6,
        49, 40.2, 53.7, 35.8,
        17.9, 17.6, 19, 36.6,
        24, 17.2, 26.5, 27.5,
        13.2, 14.7, 60.6, 10.4
    ) * 10^6
    vec_bin_count <- ceiling(vec_chromosome_bp / CN_bin_length)
    vec_centromere_location <- round(vec_centromere_bp / CN_bin_length)
    #   Set up table of chromosome bin counts and centromere locations
    TABLE_CHROMOSOME_CN_INFO <- data.frame(vec_chromosome_name, vec_bin_count, vec_centromere_location)
    columns <- c("Chromosome", "Bin_count", "Centromere_location")
    colnames(TABLE_CHROMOSOME_CN_INFO) <- columns
    #---------------------Build model input file for population dynamics
    vec_time_points <- table_population_dynamics[, 1]
    vec_cell_count <- table_population_dynamics[, 2]
    if (age_birth_unit == "day") {
        vec_time_points <- 1 * vec_time_points
    } else {
        if (age_birth_unit == "week") {
            vec_time_points <- 7 * vec_time_points
        } else {
            if (age_birth_unit == "month") {
                vec_time_points <- 30 * vec_time_points
            } else {
                if (age_birth_unit == "year") {
                    vec_time_points <- 365 * vec_time_points
                }
            }
        }
    }
    columns <- c("Age_in_day", "Total_cell_count")
    TABLE_POPULATION_DYNAMICS <- data.frame(vec_time_points, vec_cell_count)
    colnames(TABLE_POPULATION_DYNAMICS) <- columns
    #--------------------Build model input file for sampling information
    TABLE_SAMPLING_INFO <- Table_sample
    if (age_birth_unit == "day") {
        TABLE_SAMPLING_INFO$T_sample <- 1 * TABLE_SAMPLING_INFO$Age_sample
    } else {
        if (age_birth_unit == "week") {
            TABLE_SAMPLING_INFO$T_sample <- 7 * TABLE_SAMPLING_INFO$Age_sample
        } else {
            if (age_birth_unit == "month") {
                TABLE_SAMPLING_INFO$T_sample <- 30 * TABLE_SAMPLING_INFO$Age_sample
            } else {
                if (age_birth_unit == "year") {
                    TABLE_SAMPLING_INFO$T_sample <- 365 * TABLE_SAMPLING_INFO$Age_sample
                }
            }
        }
    }
    #------------------------------------Output the model variable files
    MODEL_VARIABLES <- list()
    MODEL_VARIABLES$general_variables <- TABLE_VARIABLES
    MODEL_VARIABLES$selection_model <- TABLE_SELECTION
    MODEL_VARIABLES$cn_info <- TABLE_CHROMOSOME_CN_INFO
    MODEL_VARIABLES$population_dynamics <- TABLE_POPULATION_DYNAMICS
    MODEL_VARIABLES$sampling_info <- TABLE_SAMPLING_INFO
    MODEL_VARIABLES$gc_and_mappability <- gc
    return(MODEL_VARIABLES)
}
