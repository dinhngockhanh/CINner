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
    global driver_library


    global initial_ploidy_chrom initial_ploidy_allele initial_ploidy_block
    global initial_driver_count initial_driver_map initial_DNA_length initial_selection_rate initial_prob_new_drivers
    global initial_clonal_ID initial_population initial_N_clones
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
%---Input table of CN profiles for the initial population
    filename                                = [model '-input-initial-cn-profiles.csv'];
    TABLE_INITIAL_COPY_NUMBER_PROFILES      = readtable(filename,'Delimiter',',');
%---Input other information for the initial population
    filename                                = [model '-input-initial-others.csv'];
    TABLE_INITIAL_OTHERS                    = readtable(filename,'Delimiter',',');
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
    if isempty(TABLE_CANCER_GENES)
        driver_library                      = {};
    else
        driver_library                      = TABLE_CANCER_GENES;
    end
%-----------------------------------Set up initial state for simulations
%   Get number of clones in the initial population
    initial_N_clones                        = size(TABLE_INITIAL_OTHERS,1);
    vec_header                              = TABLE_INITIAL_COPY_NUMBER_PROFILES.Properties.VariableNames;
%---Initialize ID and population for each clone
    initial_clonal_ID                       = [1:initial_N_clones];
    initial_population                      = zeros(1,initial_N_clones);
    for clone=1:initial_N_clones
%       Get driver profile of this clone
        loc                                 = find(TABLE_INITIAL_OTHERS.Clone==clone);
        population                          = TABLE_INITIAL_OTHERS.Cell_count(loc);
        initial_population(clone)           = population;
    end
%---Initialize the genotypes for each clone
    initial_ploidy_chrom                    = cell(1,initial_N_clones);
    initial_ploidy_allele                   = cell(1,initial_N_clones);
    initial_ploidy_block                    = cell(1,initial_N_clones);
    initial_driver_count                    = zeros(1,initial_N_clones);
    initial_driver_map                      = cell(1,initial_N_clones);
    initial_DNA_length                      = cell(1,initial_N_clones);
    initial_selection_rate                  = zeros(1,initial_N_clones);
    initial_prob_new_drivers                = zeros(1,initial_N_clones);
%---Set up the initial clones' CN genotypes
    for clone=1:initial_N_clones
%       Extract mini table for the CN genotypes of this clone
        text_clone_ID                       = ['Clone_' num2str(clone)];
        vec_loc                             = [1 2 find(contains(vec_header,text_clone_ID))];
        CLONE_INITIAL_COPY_NUMBER_PROFILES  = TABLE_INITIAL_COPY_NUMBER_PROFILES(:,vec_loc);
%       Set up clone's CN genotype
        ploidy_chrom                        = zeros(1,N_chromosomes);
        ploidy_block                        = cell(1,1);
        ploidy_allele                       = cell(1,1);
        for chrom=1:N_chromosomes
%           Get CN genotype for this chromosome
            CHROM_COPY_NUMBER_PROFILES      = CLONE_INITIAL_COPY_NUMBER_PROFILES(find(CLONE_INITIAL_COPY_NUMBER_PROFILES.Chromosome==chrom),:);
%           Clean CN genotype of unnecessary strands
            vec_delete                      = [];
            for column=3:size(CHROM_COPY_NUMBER_PROFILES,2)
                if all(strcmp('NA',table2cell(CHROM_COPY_NUMBER_PROFILES(:,column))))
                    vec_delete              = [vec_delete column];
                end
            end
            CHROM_COPY_NUMBER_PROFILES(:,vec_delete)    = [];
%           Update the strand count for each chromosome
            no_strands                      = size(CHROM_COPY_NUMBER_PROFILES,2)-2;
            ploidy_chrom(chrom)             = no_strands;
%           Update the CN count and allele info for each chrosomome strand
            for strand=1:no_strands
                no_blocks                   = vec_CN_block_no(chrom);
                strand_ploidy_block         = zeros(1,no_blocks);
                strand_ploidy_allele        = zeros(1,no_blocks);
                for block=1:no_blocks
                    row                     = find(CHROM_COPY_NUMBER_PROFILES.Bin==block);
                    col                     = strand+2;
                    vec_allele              = CHROM_COPY_NUMBER_PROFILES{row,col}{1};
                    if strcmp('NA',vec_allele)
                        strand_ploidy_block(block)              = 0;
                    else
                        strand_ploidy_block(block)              = length(vec_allele);
                        for unit=1:length(vec_allele)
                            strand_ploidy_allele(unit,block)    = double(vec_allele(unit))-64;
                        end
                    end
                end
                ploidy_block{chrom,strand}  = strand_ploidy_block;
                ploidy_allele{chrom,strand} = strand_ploidy_allele;
            end
        end
%       Store the clone's CN profiles
        initial_ploidy_chrom{clone}         = ploidy_chrom;
        initial_ploidy_allele{clone}        = ploidy_allele;
        initial_ploidy_block{clone}         = ploidy_block;
    end
%---Set up the initial clones' driver profiles
    for clone=1:initial_N_clones
%       Get driver profile of this clone
        loc                                 = find(TABLE_INITIAL_OTHERS.Clone==clone);
        all_drivers                         = TABLE_INITIAL_OTHERS.Drivers{loc};
        list_drivers                        = split(all_drivers,';');
%       Update the driver count for this clone
        initial_driver_count(clone)   = length(list_drivers);
%       Update the driver map for this clone
        driver_map                          = [];
        for driver=1:length(list_drivers)
            driver_info                     = split(list_drivers{driver},'_');
%           Get driver's ID, strand and unit
            driver_ID                       = driver_info{1};
            driver_strand                   = str2num(extractAfter(driver_info{2},'strand'));
            driver_unit                     = str2num(extractAfter(driver_info{3},'unit'));
%           Get driver's chromosome and block
            driver_loc                      = find(strcmp(driver_library.Gene_ID,driver_ID));
            driver_chrom                    = driver_library.Chromosome(driver_loc);
            driver_block                    = driver_library.Bin(driver_loc);
            driver_map(end+1,:)             = [driver_loc driver_chrom driver_strand driver_block driver_unit];
        end
        initial_driver_map{clone}           = driver_map;
    end
%---Set up the initial clones' DNA length
    for clone=1:initial_N_clones
        ploidy_chrom                        = initial_ploidy_chrom{clone};
        ploidy_block                        = initial_ploidy_block{clone};
        DNA_length                          = 0;
        for chrom=1:N_chromosomes
            for strand=1:ploidy_chrom(chrom)
                DNA_length                  = DNA_length+sum(ploidy_block{chrom,strand});
            end
        end
        DNA_length                          = size_CN_block_DNA*DNA_length;
        prob_new_drivers                    = 1-poisspdf(0,rate_driver*DNA_length);
        initial_DNA_length{clone}           = DNA_length;
        initial_prob_new_drivers(clone)     = prob_new_drivers;
    end
%---Set up the initial clones' selection rates
    for clone=1:initial_N_clones
        ploidy_chrom                        = initial_ploidy_chrom{clone};
        ploidy_block                        = initial_ploidy_block{clone};
        driver_count                        = initial_driver_count(clone);
        driver_map                          = initial_driver_map{clone};
        selection_rate                      = SIMULATOR_FULL_PHASE_1_selection_rate(driver_count,driver_map,ploidy_chrom,ploidy_block);
        initial_selection_rate(clone)       = selection_rate;
    end
%---Set up total population dynamics as function of age (in days)
    vec_age_in_days                         = 365*(TABLE_POPULATION_DYNAMICS.Age_in_year)';
    vec_total_cell_count                    = (TABLE_POPULATION_DYNAMICS.Total_cell_count)';
    func_expected_population                = @(time) interp1([-100*365 vec_age_in_days 200*365],[vec_total_cell_count(1) vec_total_cell_count vec_total_cell_count(end)],time);
%---Set up event rate as function of time (in days)
    cell_lifespan                           = TABLE_VARIABLES.Value(find(strcmp(TABLE_VARIABLES.Variable,'cell_lifespan')));
    func_event_rate                         = @(time) 1/cell_lifespan;
end
