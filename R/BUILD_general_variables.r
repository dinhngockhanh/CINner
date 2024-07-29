#' Produce simulations and outputs files as requested by user
#'
#' 
#' @description
#' `BUILD_general_variables` returns an object containing information about all the variables needed to produce simulations and output files specificed by the user. 
#' 
#' 
#' @param model_name A string representing the name to be given to the model. Default is 'MODEL'.
#' @param cell_lifespan An integer number representing the mean lifespan of a cell. Each cell's lifespan is exponentially distributed. Default is 4. 
#' @param cell_prob_division A float representing the probability of division of a cell. Default is NA. If NA, it overrides the selection model and the population follows neutral exponential growth.
#' @param T_0 A list consisting of an integer followed by a string specifying the start of the simulation and the unit of time eg. `list(0, 'days')`. Default is `list(0, "year")`.
#' @param T_end A list consisting of an integer followed by a string specifying the end of the simulation and the unit of time eg. `list(100, 'days')`. Default is `list(100, "year")`.
#' @param T_tau_step A numerical value represeting the time step for the tau-leaping algorithm for the simulation. Default is Inf.
#' @param Table_sample A dataframe specifying the information about sampled cells. Must have 'Age_sample' column, which is the age at which the sample is taken, usually the `T_end`. Default is `data.frame()`. 
#' @param table_population_dynamics A matrix consisting of the start and end time of the simulation as well as the initial and final population size??? Default is `matrix(0, ncol = 2, nrow = 2)`.
#' @param Max_cell_count An integer value specifying the maximum number of cells in the simulation??? Default is Inf.
#' @param Min_cell_count An integer value specifying the minimum number of cells in the simulation??? Default is 0.
#' @param Max_events An integer value specifying the maximum number of mutation events??? Default is Inf.
#' @param CN_bin_length An integer value specifying the number of basepairs belongs in a single bin of a chomosome. Default is 500000.
#' @param CN_arm_level ??? Default is FALSE.
#' @param alpha_aneuploidy An integer value representing the rate of chromosomal missegregations in WGD cells as compared to non-WGD cells. Default is 1.
#' @param prob_CN_whole_genome_duplication A numerical value representing the probability for a cell division to harbor a WGD event. Default is 0.
#' @param prob_CN_missegregation A numerical value representing the probability of whole chromosome missegregation. Default is 0.
#' @param prob_CN_chrom_arm_missegregation A numerical value representing the probability of chromosome arm missegregation. Default is 0.
#' @param prob_CN_focal_amplification A numerical value representing the probability for a cell division to harbor a focal amplification event. Default is 0.
#' @param prob_CN_focal_deletion A numerical value representing the probability for a cell division to harbor a focal deletion event. Default is 0.
#' @param prob_CN_cnloh_interstitial A numerical value representing the probability for a cell division to harbor an interstitial CN-LOH event. Default is 0.
#' @param prob_CN_cnloh_terminal A numerical value representing the probability for a cell division to harbor a terminal CN-LOH event. Default is 0.
#' @param formula_CN_whole_genome_duplication A formula of probability for a cell division to harbor a WGD event. Default is "per_division:prob_CN_whole_genome_duplication".
#' @param formula_CN_missegregation A formula of probability for a cell division to harbor a whole chromosome missegregation event. Default is "per_division:prob_CN_missegregation".
#' @param formula_CN_chrom_arm_missegregation A formula of probability for a cell division to harbor a chromosome arm missegregation event. Default is "per_division:prob_CN_chrom_arm_missegregation".
#' @param formula_CN_focal_amplification A formula of probability for a cell division to harbor a focal amplification event. Default is "per_division:prob_CN_focal_amplification".
#' @param formula_CN_focal_deletion A formula of probability for a cell division to harbor a whole focal deletion event. Default is "per_division:prob_CN_focal_deletion".
#' @param formula_CN_cnloh_interstitial A formula of probability for a cell division to harbor an interstitial CN-LOH event. Default is "per_division:prob_CN_cnloh_interstitial".
#' @param formula_CN_cnloh_terminal A formula of probability for a cell division to harbor a terminal CN-LOH event. Default is "per_division:prob_CN_cnloh_terminal".
#' @param prob_neutral_CN_whole_genome_duplication A numerical value representing the probability for a cell division to harbor a neutral WGD event. Default is 0.
#' @param prob_neutral_CN_missegregation A numerical value representing the probability for a cell division to harbor a neutral whole chromosome missegregation event. Default is 0.
#' @param prob_neutral_CN_chrom_arm_missegregation A numerical value representing the probability for a cell division to harbor a neutral chromosome arm missegregation event. Default is 0.
#' @param prob_neutral_CN_focal_amplification A numerical value representing the probability for a cell division to harbor a neutral focal amplification event. Default is 0.
#' @param prob_neutral_CN_focal_deletion A numerical value representing the probability for a cell division to harbor a neutral focal deletion event. Default is 0.
#' @param prob_neutral_CN_cnloh_interstitial A numerical value representing the probability for a cell division to harbor a neutral interstitial CN-LOH event. Default is 0.
#' @param prob_neutral_CN_cnloh_terminal A numerical value representing the probability for a cell division to harbor a neutral terminal CN-LOH event. Default is 0.
#' @param model_CN_focal_amplification_length A string value specifying the choice of distribution for the block length of a focal amplification event. Options are 'geom' and 'beta'. Default is "geom".
#' @param prob_CN_focal_amplification_length A numerical value for the geometric parameter for the block length of a focal amplification event. Only used if a geometric distribution is specified. Default is 0.1.
#' @param prob_CN_focal_amplification_length_shape_1 A numerical value for the beta parameter shape 1 for the block length of a focal amplification event. Only used if a beta distribution is specified. Default is 0.
#' @param prob_CN_focal_amplification_length_shape_2 A numerical value for the beta parameter shape 2 for the block length of a focal amplification event. Only used if a beta distribution is specified. Default is 0.
#' @param model_CN_focal_deletion_length A string value specifying the choice of distribution for the block length of a focal deletion event. Options are 'geom' and 'beta'. Default is "geom".
#' @param prob_CN_focal_deletion_length_shape_1 A numerical value for the beta parameter shape 1 for the block length of a focal deletion event. Only used if a beta distribution is specified. Default is 0.
#' @param prob_CN_focal_deletion_length_shape_2 A numerical value for the beta parameter shape 2 for the block length of a focal deletion event. Only used if a beta distribution is specified. Default is 0.
#' @param prob_CN_focal_deletion_length A numerical value for the geometric parameter for the block length of a focal deletion event. Only used if a geometric distribution is specified. Default is 0.1.
#' @param model_CN_cnloh_interstitial_length A string value specifying the choice of distribution for the block length of an interstitial CN-LOH event. Options are 'geom' and 'beta'. Default is "geom".
#' @param prob_CN_cnloh_interstitial_length A numerical value for the geometric parameter for the block length of an interstitial CN-LOH deletion event. Only used if a geometric distribution is specified. Default is 0.1.
#' @param prob_CN_cnloh_interstitial_length_shape_1 A numerical value for the beta parameter shape 1 for the block length of an interstitial CN-LOH event. Only used if a beta distribution is specified. Default is 0.
#' @param prob_CN_cnloh_interstitial_length_shape_2 A numerical value for the beta parameter shape 2 for the block length of an interstitial CN-LOH event. Only used if a beta distribution is specified. Default is 0.
#' @param model_CN_cnloh_terminal_length A string value specifying the choice of distribution for the block length of a terminal CN-LOH event. Options are 'geom' and 'beta'. Default is "geom".
#' @param prob_CN_cnloh_terminal_length A numerical value for the geometric parameter for the block length of a terminal CN-LOH deletion event. Only used if a geometric distribution is specified. Default is 0.1.
#' @param prob_CN_cnloh_terminal_length_shape_1 A numerical value for the beta parameter shape 1 for the block length of a terminal CN-LOH event. Only used if a beta distribution is specified. Default is 0.
#' @param prob_CN_cnloh_terminal_length_shape_2 A numerical value for the beta parameter shape 2 for the block length of a terminal CN-LOH event. Only used if a beta distribution is specified. Default is 0.
#' @param rate_driver A numerical value for the poisson rate of getting new driver mutations. Default is 0.
#' @param rate_passenger A numerical value for the poisson rate of getting a new passenger mutation. Default is 0.
#' @param selection_model A string specifying the selection model to be used. ???
#' @param bound_driver An integer number representing the maximum driver count in viable cells (cells exceeding this will die). Default is Inf.
#' @param bound_average_ploidy An integer number representing the maximum average ploidy across genome (cells exceeding this will die). Default is Inf.
#' @param bound_homozygosity An integer number representing the maximum number of bins under homozygosity (cells exceeding this will die). Default is 0.
#' @param SFS_totalsteps An integer number representing the bin count in SFS data. Default is 25.
#' @param prob_coverage A numerical value representing the mean coverage depth. Default is 0.05.
#' @param alpha_coverage A numerical value representing the alpha parameter for coverage depth. Default is 0.7.
#' @param lower_limit_cell_counts A numerical value representing the lower limit of cell counts for mutations to be detected. Default is 0.
#' @param lower_limit_alt_counts A numerical value representing the lower limit of alternate read counts for mutations to be detected. Default is 3.
#' @param lower_limit_tot_counts A numerical value representing the lower limit of total read counts for mutations to be detected. Default is 0.
#' @param gc ???
#' @param gc_slope A numerical value for the slope of the linear GC model. Default is 1.2.
#' @param gc_int A numerical value for the intercept of the linear GC model. Default is 0.
#' @param sigma1 A numerical value for the gamma scale for read depth noise. Default is 0.1.
#' @param num_reads A numerical value for the number of reads per cell. Default is 1e6. 
#' 
#' @examples
#' 
#' cell_lifespan <- 1
#' T_0 <- list(0, "day")
#' T_end <- list(300, "day")
#' Table_sample <- data.frame(Sample_ID = c("SA"), Cell_count = c(Inf), Age_sample = c(T_end[[1]]))
#' T_tau_step <- cell_lifespan / 2
#' CN_bin_length <- 500000
#' 
#' selection_model <- "chrom-arm-selection"
#' 
#' prob_CN_whole_genome_duplication <- 0
#' prob_CN_missegregation <- 0
#' prob_CN_chrom_arm_missegregation <- 0
#' prob_CN_focal_amplification <- 0
#' prob_CN_focal_deletion <- 0
#' prob_CN_cnloh_interstitial <- 0
#' prob_CN_cnloh_terminal <- 0
#' model_CN_focal_amplification_length <- "beta"
#' model_CN_focal_deletion_length <- "beta"
#' prob_CN_focal_amplification_length_shape_1 <- 0.758304780825031
#' prob_CN_focal_amplification_length_shape_2 <- 5.33873409782625
#' prob_CN_focal_deletion_length_shape_1 <- 0.814054548726361
#' prob_CN_focal_deletion_length_shape_2 <- 6.16614890284825
#' prob_CN_cnloh_interstitial_length <- 0.005
#' prob_CN_cnloh_terminal_length <- 0.005
#' rate_driver <- 0
#' rate_passenger <- 1e-11
#' 
#' bound_driver <- 3
#' bound_average_ploidy <- 6
#' bound_maximum_CN <- 8
#' bound_homozygosity <- 0
#' 
#' vec_time <- T_0[[1]]:T_end[[1]]
#' vec_cell_count <- rep(1000,length(vec_time))
#' table_population_dynamics <- cbind(vec_time, vec_cell_count)
#'
#' gc <- read.csv(file = system.file("extdata", "gc_map_500kb.csv", package = "CINner"))
#' gc_slope <- 1.2
#' gc_int <- 0
#' sigma1 <- 0.02642392
#' num_reads <- 3906632
#'
#' model_variables_base <- BUILD_general_variables(
#'    cell_lifespan = cell_lifespan,
#'    T_0 = T_0, T_end = T_end, T_tau_step = T_tau_step,
#'    Table_sample = Table_sample,
#'    CN_bin_length = CN_bin_length,
#'    prob_CN_whole_genome_duplication = prob_CN_whole_genome_duplication,
#'    prob_CN_missegregation = prob_CN_missegregation,
#'    prob_CN_chrom_arm_missegregation = prob_CN_chrom_arm_missegregation,
#'    prob_CN_focal_amplification = prob_CN_focal_amplification,
#'    prob_CN_focal_deletion = prob_CN_focal_deletion,
#'    prob_CN_cnloh_interstitial = prob_CN_cnloh_interstitial,
#'    prob_CN_cnloh_terminal = prob_CN_cnloh_terminal,
#'    model_CN_focal_amplification_length = model_CN_focal_amplification_length,
#'    model_CN_focal_deletion_length = model_CN_focal_deletion_length,
#'    prob_CN_focal_amplification_length_shape_1 = prob_CN_focal_amplification_length_shape_1,
#'    prob_CN_focal_amplification_length_shape_2 = prob_CN_focal_amplification_length_shape_2,
#'    prob_CN_focal_deletion_length_shape_1 = prob_CN_focal_deletion_length_shape_1,
#'    prob_CN_focal_deletion_length_shape_2 = prob_CN_focal_deletion_length_shape_2,
#'    prob_CN_cnloh_interstitial_length = prob_CN_cnloh_interstitial_length,
#'    prob_CN_cnloh_terminal_length = prob_CN_cnloh_terminal_length,
#'    rate_driver = rate_driver,
#'    rate_passenger = rate_passenger,
#'    selection_model = selection_model,
#'    bound_driver = bound_driver,
#'    bound_average_ploidy = bound_average_ploidy,
#'    bound_maximum_CN = bound_maximum_CN,
#'    bound_homozygosity = bound_homozygosity,
#'    table_population_dynamics = table_population_dynamics,
#'    gc = gc,
#'    gc_slope = gc_slope,
#'    gc_int = gc_int,
#'    sigma1 = sigma1,
#'    num_reads = num_reads)
#' 
#' @export
BUILD_general_variables <- function(model_name = "MODEL",
                                    cell_lifespan = 4,
                                    cell_prob_division = NA,
                                    T_0 = list(0, "year"),
                                    T_end = list(100, "year"),
                                    T_tau_step = Inf,
                                    Table_sample = data.frame(),
                                    table_population_dynamics = matrix(0, ncol = 2, nrow = 2),
                                    Max_cell_count = Inf,
                                    Min_cell_count = 0,
                                    Max_events = Inf,
                                    CN_bin_length = 500000,
                                    CN_arm_level = FALSE,
                                    # model_CN_whole_genome_duplication = "per_division",
                                    # model_CN_missegregation = "per_division",
                                    # model_CN_chrom_arm_missegregation = "per_division",
                                    # model_CN_focal_amplification = "per_division",
                                    # model_CN_focal_deletion = "per_division",
                                    # model_CN_cnloh_interstitial = "per_division",
                                    # model_CN_cnloh_terminal = "per_division",
                                    alpha_aneuploidy = 1,
                                    prob_CN_whole_genome_duplication = 0,
                                    prob_CN_missegregation = 0,
                                    prob_CN_chrom_arm_missegregation = 0,
                                    prob_CN_focal_amplification = 0,
                                    prob_CN_focal_deletion = 0,
                                    prob_CN_cnloh_interstitial = 0,
                                    prob_CN_cnloh_terminal = 0,
                                    formula_CN_whole_genome_duplication = "per_division:prob_CN_whole_genome_duplication",
                                    formula_CN_missegregation = "per_division:prob_CN_missegregation",
                                    formula_CN_chrom_arm_missegregation = "per_division:prob_CN_chrom_arm_missegregation",
                                    formula_CN_focal_amplification = "per_division:prob_CN_focal_amplification",
                                    formula_CN_focal_deletion = "per_division:prob_CN_focal_deletion",
                                    formula_CN_cnloh_interstitial = "per_division:prob_CN_cnloh_interstitial",
                                    formula_CN_cnloh_terminal = "per_division:prob_CN_cnloh_terminal",
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
                                    bound_driver = Inf,
                                    bound_average_ploidy = Inf,
                                    bound_homozygosity = 0,
                                    bound_maximum_CN = Inf,
                                    bound_maximum_CN_normalized = Inf,
                                    bound_WGD = Inf,
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
    TABLE_VARIABLES[N_row, ] <- c("Max_cell_count", Max_cell_count, "cell count", "Maximum cell count, simulation is ended if violated (Inf if no condition)")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("Min_cell_count", Min_cell_count, "cell count", "Minimum cell count, simulation is rejected if violated (0 if no condition)")
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("Max_events", Max_events, "event count", "Maximum event count, simulation is ended if violated (Inf if no condition)")
    #   Set up constant cell lifespan
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("cell_lifespan", cell_lifespan, "day", "Lifespan of one cell")
    #   Set up constant cell division probability
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("cell_prob_division", cell_prob_division, "", "Division probability of one cell (if NA: overrides the selection model, population follows neutral exponential growth)")
    #   Set up CN bin width
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("size_CN_block_DNA", CN_bin_length, "bp", "CN bin width")

    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("alpha_aneuploidy", alpha_aneuploidy, "", "Degree of aneuploidy associated with increased CNA rates depending on genotype")

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
    TABLE_SELECTION[N_row, ] <- c("bound_WGD", bound_WGD, "", "Maximum WGD allowed for a cell (cells are not allowed to exceed this)")
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
    if (!CN_arm_level) {
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
    } else {
        vec_bin_count <- rep(2, length(vec_chromosome_name))
        vec_centromere_location <- rep(1, length(vec_chromosome_name))
    }
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
