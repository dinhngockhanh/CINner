%=============================SET ALL PARAMETERS REQUIRED FOR SIMULATION
function SIMULATOR_VARIABLES_for_simulation(model)
    global func_expected_population func_event_rate
    global T_start_time T_end_time age_birth age_end Population_end Max_events
    global SFS_totalsteps N_sample
    global bound_ploidy


    global prob_CN_whole_genome_duplication
    global prob_CN_missegregation
    global prob_CN_chrom_arm_missegregation
    global prob_CN_focal_amplification prob_CN_focal_amplification_length
    global prob_CN_focal_deletion prob_CN_focal_deletion_length
    global prob_CN_cnloh_interstitial prob_CN_cnloh_interstitial_length
    global prob_CN_cnloh_terminal prob_CN_cnloh_terminal_length


    global bound_driver rate_driver rate_passenger
    global growth_model carrying_capacity rate_base_lifetime rate_selection_recurrently_mutated_genes rate_selection_recurrently_CNA_genes prob_division
    global logistic_initial logistic_slope logistic_halftime
    global level_purity prob_coverage alpha_coverage lower_limit_cell_counts lower_limit_alt_counts lower_limit_tot_counts
    global vec_SEER_age
    global T_tau_step
    global N_chromosomes size_CN_block_DNA vec_CN_block_no vec_centromere_location
    global driver_library driver_order
%======================================================GENERAL VARIABLES
%   T_end_time                          : Final time of simulation
%   Population_end                      : Target for final population size
%   Max_events                          : Limit on number of events
%   N_sample                            : Number of cells in sequencing sample
%   SFS_totalsteps                      : Number of steps to divide SFS frequencies in [0,1]
%======================================CHROMOSOME PLOIDY AND COPY NUMBER
%   bound_ploidy                        : Bound on ploidy in one cell
%   prob_CN_whole_genome_duplication    : Probability of whole genome amplification
%   prob_CN_missegregation              : Probability of chromosome missegregation
%   prob_CN_chrom_arm_missegregation    : Probability of chromosome arm missegregation
%   prob_CN_focal_amplification         : Probability of focal CN event in a chromosome
%   prob_CN_focal_amplification_length      : Geometric parameter for the block length of focal CN event
%=========================================DRIVER AND PASSENGER MUTATIONS
%   bound_driver                        : Bound on driver count in one cell
%   rate_driver                         : Poisson rate of new driver mutations
%   rate_passenger                      : Poisson rate of new passenger mutations
%=================================================CELL LIFETIME AND FATE
%   growth_model                        : Model of growth
%   func_expected_population            : Function for expected population size, given the current time point
%   carrying_capacity                   : Carrying capacity
%   logistic_initial                    : Initial division rate for logistic growth model
%   logistic_slope                      : Slope for logistic growth model
%   logistic_halftime                   : Halftime for logistic growth model
%   func_event_rate                     : Function for base rate for a cell's lifetime, given the current time point
%   rate_selection_recurrently_mutated_genes    : Selection rate for a cell's life time with recurrently mutated drivers
%   rate_selection_recurrently_CNA_genes        :Selection rate for a cell's life time with recurrently CNA drivers
%========================================SAMPLING AND SEQUENCING EFFECTS
%   level_purity                        : Assumed level of purity in cancer samples
%   prob_coverage                       : Sequencing coverage (probability that a read is recorded)
%   alpha_coverage                      : Shape parameter of sequencing coverage
%   lower_limit_cell_counts             : Lower limit on alternate cell counts for SFS
%   lower_limit_alt_counts              : Lower limit on alternate read counts for SFS
%   lower_limit_tot_counts              : Lower limit on total read counts for SFS
%--------------------------------------Set up variables for cancer model
%---Input table of variables from file
    filename                                = [model '-input-variables.csv'];
    TABLE_VARIABLES                         = readtable(filename,'Delimiter',',');
%---Input table of total population dynamics
    filename                                = [model '-input-population-dynamics.csv'];
    TABLE_POPULATION_DYNAMICS               = readtable(filename,'Delimiter',',');
%---Input table of chromosome bin counts and centromere locations
    filename                                = [model '-input-copy-number-blocks.csv'];
    TABLE_CHROMOSOME_CN_INFO                = readtable(filename,'Delimiter',',');
%---Input table of mutational and CNA genes
    filename                                = [model '-input-cancer-genes.csv'];
    TABLE_CANCER_GENES                      = readtable(filename,'Delimiter',',');





%---Set up individual variables from table
    age_birth                               = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'age_birth')));
    age_end                                 = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'age_end')));
    T_start_time                            = age_birth*365;
    T_end_time                              = age_end*365;
    Population_end                          = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'Population_end')));
    Max_events                              = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'Max_events')));
    N_sample                                = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'N_sample')));
    SFS_totalsteps                          = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'SFS_totalsteps')));
    rate_driver                             = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'rate_driver')));
    bound_driver                            = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'bound_driver')));
    bound_ploidy                            = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'bound_ploidy')));

    prob_CN_whole_genome_duplication        = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_whole_genome_duplication')));
    prob_CN_missegregation                  = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_missegregation')));
    prob_CN_chrom_arm_missegregation        = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_chrom_arm_missegregation')));
    prob_CN_focal_amplification             = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_focal_amplification')));
    prob_CN_focal_deletion                  = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_focal_deletion')));
    prob_CN_cnloh_interstitial              = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_cnloh_interstitial')));
    prob_CN_cnloh_terminal                  = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_cnloh_terminal')));

    prob_CN_focal_amplification_length      = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_focal_amplification_length')));
    prob_CN_focal_deletion_length           = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_focal_deletion_length')));
    prob_CN_cnloh_interstitial_length       = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_cnloh_interstitial_length')));
    prob_CN_cnloh_terminal_length           = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_CN_cnloh_terminal_length')));

    rate_passenger                          = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'rate_passenger')));
    prob_coverage                           = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'prob_coverage')));
    alpha_coverage                          = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'alpha_coverage')));
    lower_limit_cell_counts                 = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'lower_limit_cell_counts')));
    lower_limit_alt_counts                  = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'lower_limit_alt_counts')));
    lower_limit_tot_counts                  = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'lower_limit_tot_counts')));
    T_tau_step                              = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'T_tau_step')));



%   level_purity is only for SIMULATOR_CLONAL and SIMULATOR_ODE
    level_purity                            = 1;



%---Set up variables for copy number model
    N_chromosomes                           = size(TABLE_CHROMOSOME_CN_INFO,1);
    size_CN_block_DNA                       = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'size_CN_block_DNA')));
    vec_CN_block_no                         = TABLE_CHROMOSOME_CN_INFO.Bin_count';
    vec_centromere_location                 = TABLE_CHROMOSOME_CN_INFO.Centromere_location';
%---Set up mutational and CNA driver library (without selection rates)
    driver_order                            = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'driver_order')));
    if isnan(driver_order)
        driver_order                        = {};
    else
%       ...
%       ...
%       ...
%       ...
%       ...
%       ...
%       ...
%       ...
%       ...
%       ...
    end
    if isempty(TABLE_CANCER_GENES)
        driver_library                      = {};
    else
        driver_library                      = {};
        for row=1:size(TABLE_CANCER_GENES)
            driver_library{row,1}           = TABLE_CANCER_GENES.Gene_ID{row};
            driver_library{row,2}           = TABLE_CANCER_GENES.Gene_role{row};
            driver_library{row,3}           = TABLE_CANCER_GENES.Chromosome(row);
            driver_library{row,4}           = TABLE_CANCER_GENES.Bin(row);
            driver_library{row,5}           = TABLE_CANCER_GENES.Affected_donor_count(row);
        end
    end
%---Set up total population dynamics as function of age (in days)
    vec_age_in_days                         = 365*(TABLE_POPULATION_DYNAMICS.Age_in_year)';
    vec_total_cell_count                    = (TABLE_POPULATION_DYNAMICS.Total_cell_count)';
    func_expected_population                = @(time) interp1([-100*365 vec_age_in_days 200*365],[vec_total_cell_count(1) vec_total_cell_count vec_total_cell_count(end)],time);
%---Set up event rate as function of time (in days)
    cell_lifespan                           = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'cell_lifespan')));
    func_event_rate                         = @(time) 1/cell_lifespan;
end
