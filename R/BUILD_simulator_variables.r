BUILD_simulator_variables <- function(model_name                        = 'MODEL',
                                      cell_lifespan                     = 4,
                                      T_0                               = list(0,'year'),
                                      T_end                             = list(100,'year'),
                                      Population_end                    = Inf,
                                      Max_events                        = Inf,
                                      CN_bin_length                     = 500000,
                                      prob_CN_whole_genome_duplication  = 0,
                                      prob_CN_missegregation            = 0,
                                      prob_CN_chrom_arm_missegregation  = 0,
                                      prob_CN_focal_amplification       = 0,
                                      prob_CN_focal_deletion            = 0,
                                      prob_CN_cnloh_interstitial        = 0,
                                      prob_CN_cnloh_terminal            = 0,
                                      prob_CN_focal_amplification_length= 0.1,
                                      prob_CN_focal_deletion_length     = 0.1,
                                      prob_CN_cnloh_interstitial_length = 0.1,
                                      prob_CN_cnloh_terminal_length     = 0.1){
#---------------------------Build model input file for general variables
    columns                             <- c('Variable','Value','Unit','Note')
    TABLE_VARIABLES                     <- data.frame(matrix(nrow=0,ncol=length(columns)))
    colnames(TABLE_VARIABLES)           <- columns
    N_row                               <- 0
#   Set up the start time of simulations
    age_birth                           <- T_0[[1]]
    age_birth_unit                      <- T_0[[2]]
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('age_birth',age_birth,age_birth_unit,'Age when simulation starts')
    if (age_birth_unit=='day'){
        T_start_time                    <- age_birth
    }else{if (age_birth_unit=='week'){
        T_start_time                    <- 7*age_birth
    }else{if (age_birth_unit=='month'){
        T_start_time                    <- 30*age_birth
    }else{if (age_birth_unit=='year'){
        T_start_time                    <- 365*age_birth
    }}}}
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('T_start_time',T_start_time,'day','Age when simulation starts (for internal use)')
#   Set up the end time of simulations
    age_end                             <- T_end[[1]]
    age_end_unit                        <- T_end[[2]]
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('age_end',age_end,age_end_unit,'Age when simulation stops')
    if (age_end_unit=='day'){
        T_end_time                      <- age_end
    }else{if (age_end_unit=='week'){
        T_end_time                      <- 7*age_end
    }else{if (age_end_unit=='month'){
        T_end_time                      <- 30*age_end
    }else{if (age_end_unit=='year'){
        T_end_time                      <- 365*age_end
    }}}}
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('T_end_time',T_end_time,'day','Age when simulation stops (for internal use)')
#   Set up condition to end simulation prematurely based on population size or event count
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('Population_end',Population_end,'cell count','Condition for ending simulation (Inf if no condition)')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('Max_events',Max_events,'event count','Condition for ending simulation (Inf if no condition)')
#   Set up cell lifespan
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('cell_lifespan',cell_lifespan,'day','Lifespan of one cell')
#   Set up CN bin width
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('size_CN_block_DNA',CN_bin_length,'bp','CN bin width')
#   Set up CN event probabilities
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_whole_genome_duplication',prob_CN_whole_genome_duplication,'per cell division','Probability for a cell division to harbor a WGD event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_missegregation',prob_CN_missegregation,'per cell division','Probability for a cell division to harbor a chromosome mis-segregation event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_chrom_arm_missegregation',prob_CN_chrom_arm_missegregation,'per cell division','Probability for a cell division to harbor a chromosome-arm mis-segregation event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_amplification',prob_CN_focal_amplification,'per cell division','Probability for a cell division to harbor a focal amplification event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_deletion',prob_CN_focal_deletion,'per cell division','Probability for a cell division to harbor a focal deletion event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_interstitial',prob_CN_cnloh_interstitial,'per cell division','Probability for a cell division to harbor an interstitial CN-LOH event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_terminal',prob_CN_cnloh_terminal,'per cell division','Probability for a cell division to harbor a terminal CN-LOH event')
#   Set up geometric parameters for lengths of local CN events
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_amplification_length',prob_CN_focal_amplification_length,,'Geometric parameter for the block length of a focal amplification event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_deletion_length',prob_CN_focal_deletion_length,,'Geometric parameter for the block length of a focal deletion event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_interstitial_length',prob_CN_cnloh_interstitial_length,,'Geometric parameter for the block length of an interstitial CN-LOH event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_terminal_length',prob_CN_cnloh_terminal_length,,'Geometric parameter for the block length of a terminal CN-LOH event')





#-----------------------Build model input file for chromosome bin counts
#-----------------------------------------------and centromere locations
    vec_chromosome_name                 <- c(1:22,'X','Y')
    vec_chromosome_bp                   <- c(248956422,     242193529,      198295559,      190214555,
                                             181538259,     170805979,      159345973,      145138636,
                                             138394717,     133797422,      135086622,      133275309,
                                             114364328,     107043718,      101991189,      90338345,
                                             83257441,      80373285,       58617616,       64444167,
                                             46709983,      50818468,       156040895,      57227415)
    vec_centromere_bp                   <- c(125,           93.3,           91,             50.4,
                                             48.4,          61,             59.9,           45.6,
                                             49,            40.2,           53.7,           35.8,
                                             17.9,          17.6,           19,             36.6,
                                             24,            17.2,           26.5,           27.5,
                                             13.2,          14.7,           60.6,           10.4)*10^6
    vec_bin_count                           <- ceiling(vec_chromosome_bp/CN_bin_length)
    vec_centromere_location             <- round(vec_centromere_bp/CN_bin_length)
#   Set up table of chromosome bin counts and centromere locations
    columns                             <- c('Chromosome','Bin_count','Centromere_location')
    TABLE_CHROMOSOME_CN_INFO            <- data.frame(vec_chromosome_name,vec_bin_count,vec_centromere_location)
    colnames(TABLE_CHROMOSOME_CN_INFO)  <- columns






    print(TABLE_VARIABLES)
    print(TABLE_CHROMOSOME_CN_INFO)
}
