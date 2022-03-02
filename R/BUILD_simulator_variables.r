
BUILD_general_variables <- function(model_name                           = 'MODEL',
                                    cell_lifespan                        = 4,
                                    T_0                                  = list(0,'year'),
                                    T_end                                = list(100,'year'),
                                    T_tau_step                           = 3,
                                    table_population_dynamics            = matrix(0,ncol=2,nrow=2),
                                    Population_end                       = Inf,
                                    Max_events                           = Inf,
                                    CN_bin_length                        = 500000,
                                    prob_CN_whole_genome_duplication     = 0,
                                    prob_CN_missegregation               = 0,
                                    prob_CN_chrom_arm_missegregation     = 0,
                                    prob_CN_focal_amplification          = 0,
                                    prob_CN_focal_deletion               = 0,
                                    prob_CN_cnloh_interstitial           = 0,
                                    prob_CN_cnloh_terminal               = 0,
                                    prob_CN_focal_amplification_length   = 0.1,
                                    prob_CN_focal_deletion_length        = 0.1,
                                    prob_CN_cnloh_interstitial_length    = 0.1,
                                    prob_CN_cnloh_terminal_length        = 0.1,
                                    rate_driver                          = 0,
                                    rate_passenger                       = 0,
                                    bound_driver                         = 3,
                                    bound_ploidy                         = 10,
                                    SFS_totalsteps                       = 25,
                                    prob_coverage                        = 0.05,
                                    alpha_coverage                       = 0.7,
                                    lower_limit_cell_counts              = 0,
                                    lower_limit_alt_counts               = 3,
                                    lower_limit_tot_counts               = 0){
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
    if (age_birth_unit!=age_end_unit){
        print('START AND END TIMES DO NOT USE THE SAME UNIT')
    }

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
#   Set up the time step for tau-leaping algorithm
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('T_tau_step',T_tau_step,'day','Time step for tau-leaping algorithm for simulation')
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
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_amplification_length',prob_CN_focal_amplification_length,'','Geometric parameter for the block length of a focal amplification event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_focal_deletion_length',prob_CN_focal_deletion_length,'','Geometric parameter for the block length of a focal deletion event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_interstitial_length',prob_CN_cnloh_interstitial_length,'','Geometric parameter for the block length of an interstitial CN-LOH event')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_CN_cnloh_terminal_length',prob_CN_cnloh_terminal_length,'','Geometric parameter for the block length of a terminal CN-LOH event')
#   Set up mutation rates
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('rate_driver',rate_driver,'per bp per cell division','Poisson rate of getting new driver mutations')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('rate_passenger',rate_passenger,'per bp per cell division','Poisson rate of getting new passenger mutations')
#   Set up upper limits for driver and local CN for viable cells
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('bound_driver',bound_driver,'driver count','Maximum driver count in viable cells (cells exceeding this will die)')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('bound_ploidy',bound_ploidy,'local CN','Maximum local CN in viable cells (cells exceeding this will die)')
#   Set up variables for sequencing
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('SFS_totalsteps',SFS_totalsteps,'','Bin count in SFS data')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('prob_coverage',prob_coverage,'','Mean coverage depth')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('alpha_coverage',alpha_coverage,'','Alpha parameter for coverage depth')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('lower_limit_cell_counts',lower_limit_cell_counts,'','Lower limit of cell counts for mutations to be detected')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('lower_limit_alt_counts',lower_limit_alt_counts,'','Lower limit of alternate read counts for mutations to be detected')
    N_row                               <- N_row+1
    TABLE_VARIABLES[N_row,]             <- c('lower_limit_tot_counts',lower_limit_tot_counts,'','Lower limit of total read counts for mutations to be detected')
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
    vec_bin_count                       <- ceiling(vec_chromosome_bp/CN_bin_length)
    vec_centromere_location             <- round(vec_centromere_bp/CN_bin_length)
#   Set up table of chromosome bin counts and centromere locations
    columns                             <- c('Chromosome','Bin_count','Centromere_location')
    TABLE_CHROMOSOME_CN_INFO            <- data.frame(vec_chromosome_name,vec_bin_count,vec_centromere_location)
    colnames(TABLE_CHROMOSOME_CN_INFO)  <- columns
#-------------------------Build model input file for population dynamics
    vec_time_points                     <- table_population_dynamics[,1]
    vec_cell_count                      <- table_population_dynamics[,2]
    if (age_birth_unit=='day'){
        vec_time_points                 <- 1*vec_time_points
    }else{if (age_birth_unit=='week'){
        vec_time_points                 <- 7*vec_time_points
    }else{if (age_birth_unit=='month'){
        vec_time_points                 <- 30*vec_time_points
    }else{if (age_birth_unit=='year'){
        vec_time_points                 <- 365*vec_time_points
    }}}}
    columns                             <- c('Age_in_day','Total_cell_count')
    TABLE_POPULATION_DYNAMICS           <- data.frame(vec_time_points,vec_cell_count)
    colnames(TABLE_POPULATION_DYNAMICS) <- columns
#----------------------------------------Output the model variable files
    MODEL_VARIABLES                     <- list()
    MODEL_VARIABLES$general_variables   <- TABLE_VARIABLES
    MODEL_VARIABLES$cn_info             <- TABLE_CHROMOSOME_CN_INFO
    MODEL_VARIABLES$population_dynamics <- TABLE_POPULATION_DYNAMICS
    return(MODEL_VARIABLES)
}

BUILD_driver_library <- function(MODEL_VARIABLES    = list(),
                                 vec_driver_genes   = c(),
                                 vec_driver_role    = c(),
                                 vec_chromosome     = -1,
                                 vec_bin            = -1,
                                 vec_driver_s       = c()){
#-------------------------------------------Input the Cancer Gene Census
    DATA_cancer_gene_census                     <- read.csv('../data/cancer_gene_census.csv')
#-----------------------------------------------Build the driver library
    columns                                     <- c('Gene_ID','Gene_role','Selective_strength')
    TABLE_CANCER_GENES                          <- data.frame(vec_driver_genes,vec_driver_role,vec_driver_s)
    colnames(TABLE_CANCER_GENES)                <- columns
    if (any(vec_chromosome==-1)){
        TABLE_CANCER_GENES$Chromosome           <- -1
    }else{
        TABLE_CANCER_GENES$Chromosome           <- vec_chromosome
    }
    if (any(vec_bin==-1)){
        TABLE_CANCER_GENES$Bin                  <- -1
    }else{
        TABLE_CANCER_GENES$Bin                  <- vec_bin
    }
#--------------------------Supplement driver library with gene locations
    if ((any(TABLE_CANCER_GENES$Chromosome==-1)==TRUE) | (any(TABLE_CANCER_GENES$Bin==-1)==TRUE)){
#       Get the CN bin length
        size_CN_block_DNA                       <- as.numeric(MODEL_VARIABLES$general_variables$Value[MODEL_VARIABLES$general_variables$Variable=='size_CN_block_DNA'])
#       Find the chromosome and bin for each driver
        for (gene in 1:length(vec_driver_genes)){
            Gene_ID                             <- vec_driver_genes[gene]
            loc                                 <- which(DATA_cancer_gene_census$Gene.Symbol==Gene_ID)
            Gene_address                        <- DATA_cancer_gene_census$Genome.Location[loc]
            Gene_chromosome                     <- as.numeric(sub(':.*','',Gene_address))
            Gene_bin                            <- round(as.numeric(sub('-.*','',sub('.*:','',Gene_address)))/size_CN_block_DNA)
            TABLE_CANCER_GENES$Chromosome[gene] <- Gene_chromosome
            TABLE_CANCER_GENES$Bin[gene]        <- Gene_bin
        }
    }
#------------------Supplement driver library with allele selection rates
#   Count the number of TSGs and ONCOGENEs
    count_TSG                                   <- sum(TABLE_CANCER_GENES$Gene_role=='TSG')
    count_ONCOGENE                              <- sum(TABLE_CANCER_GENES$Gene_role=='ONCOGENE')
#---Compute selection rates for TSGs
    list_TSG                                    <- which(TABLE_CANCER_GENES$Gene_role=='TSG')
    s_normalization                             <- 1
    for (driver in 1:length(list_TSG)){
        row                                     <- list_TSG[driver]
#       Get its selection strength
        driver_sel_rate                         <- vec_driver_s[row]
#       Compute its selection rate for WT and MUT alleles
        TABLE_CANCER_GENES$s_rate_WT[row]       <- 1/(1+driver_sel_rate)
        TABLE_CANCER_GENES$s_rate_MUT[row]      <- 1
#       Update normalizer for selection rate
        s_normalization                         <- s_normalization*(1+driver_sel_rate)
    }
#---Compute selection rates for ONCOGENEs
    s_normalization                             <- s_normalization^(1/count_ONCOGENE)
    list_ONCOGENE                               <- which(TABLE_CANCER_GENES$Gene_role=='ONCOGENE')
    for (driver in 1:length(list_ONCOGENE)){
        row                                     <- list_ONCOGENE[driver]
#       Get its selection strength
        driver_sel_rate                         <- vec_driver_s[row]
#       Compute its selection rate for WT and MUT alleles
        TABLE_CANCER_GENES$s_rate_WT[row]       <- s_normalization
        TABLE_CANCER_GENES$s_rate_MUT[row]      <- s_normalization*(1+driver_sel_rate)
    }
#----------------------------------------Output the model variable files
    MODEL_VARIABLES$driver_library              <- TABLE_CANCER_GENES
    return(MODEL_VARIABLES)
}

BUILD_initial_population <- function(MODEL_VARIABLES    = list(),
                                     cell_count         = 1,
                                     CN_arm             = cbind(c(rep('A-A',length=23),''),c(rep('B-B',length=23),'')),
                                     CN_focal           = list(),
                                     drivers            = list()){
#------------------------------Update initial clones - other information
    if (is.null(MODEL_VARIABLES$initial_others)){
        columns                                             <- c('Clone','Cell_count','Drivers')

        DRIVERS                                             <- ''
        if (length(drivers)>0){
            # for (i in 1:)
        }

        TABLE_INITIAL_OTHERS                                <- data.frame(1,cell_count,drivers)
        colnames(TABLE_INITIAL_OTHERS)                      <- columns
        I_clone                                             <- 1
    }else{
        TABLE_INITIAL_OTHERS                                <- MODEL_VARIABLES$initial_others
        I_clone                                             <- nrow(TABLE_INITIAL_OTHERS)+1
        TABLE_INITIAL_OTHERS[I_clone,]                      <- c(1,cell_count,drivers)
    }
#------------------------------------Update initial clones - CN profiles
#   Initialize table of alleles if necessary
    if (is.null(MODEL_VARIABLES$initial_cn)){
        TABLE_CHROMOSOME_CN_INFO                            <- MODEL_VARIABLES$cn_info
        vec_chrom                                           <- c()
        vec_bin                                             <- c()
        for (chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
            bin_count                                       <- TABLE_CHROMOSOME_CN_INFO$Bin_count[chrom]
            vec_chrom                                       <- c(vec_chrom,rep(TABLE_CHROMOSOME_CN_INFO$Chromosome[chrom],bin_count))
            vec_bin                                         <- c(vec_bin,(1:bin_count))
        }
        columns                                             <- c('Chromosome','Bin')
        TABLE_INITIAL_CN                                    <- data.frame(vec_chrom,vec_bin)
        colnames(TABLE_INITIAL_CN)                          <- columns
    }else{
        TABLE_INITIAL_CN                                    <- MODEL_VARIABLES$initial_cn
    }
#---Build the CN profile for the new clone - start from arm level
    max_no_strands                                          <- ncol(CN_arm)
    for (strand in 1:max_no_strands){
        TABLE_INITIAL_CN[paste('Clone_',I_clone,'_strand_',strand,sep='')]  <- 'NA'

    }
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                                               <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        centromere                                          <- TABLE_CHROMOSOME_CN_INFO$Centromere_location[i_chrom]
        for (strand in 1:max_no_strands){
            col                                             <- which(colnames(TABLE_INITIAL_CN)==paste('Clone_',I_clone,'_strand_',strand,sep=''))
            strand_allele                                   <- CN_arm[i_chrom,strand]
            if (strand_allele==''){
                next
            }
            strand_allele_left                              <- sub('-.*','',strand_allele)
            rows                                            <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin<=centromere))
            TABLE_INITIAL_CN[rows,col]                      <- strand_allele_left

            strand_allele_right                             <- sub('.*-','',strand_allele)
            rows                                            <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin>centromere))
            TABLE_INITIAL_CN[rows,col]                      <- strand_allele_right
        }
    }




#     for (i in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
#         chrom                                               <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i]
#         centromere                                          <- TABLE_CHROMOSOME_CN_INFO$Centromere_location[i]
#         for (arm in 1:2){
#             if (arm==1){
#                 vec_rows                                    <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin<=centromere))
#                 vec_alleles                                 <- CN_arm[2*i-1]
#             }else{
#                 vec_rows                                    <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin>centromere))
#                 vec_alleles                                 <- CN_arm[2*i]
#             }
#             if (nchar(vec_alleles)==0){
#                 next
#             }
#             for (strand in 1:nchar(vec_alleles)){
#                 allele                                      <- substr(vec_alleles,strand,strand)
#                 col                                         <- which(colnames(TABLE_INITIAL_CN)==paste('Clone_',I_clone,'_strand_',strand,sep=''))
#                 TABLE_INITIAL_CN[vec_rows,col]              <- allele
#             }
#         }
#     }
# #---Build the CN profile for the new clone - continue with local regions
#     if (length(CN_focal)>=1){
#         for (i in 1:length(CN_focal)){
#             chrom                                           <- CN_focal[[i]][[1]]
#             bin_start                                       <- CN_focal[[i]][[2]]
#             bin_end                                         <- CN_focal[[i]][[3]]
#             vec_alleles                                     <- CN_focal[[i]][[4]]
#             vec_rows                                        <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin>=bin_start)&(TABLE_INITIAL_CN$Bin<=bin_end))
#             if ((vec_alleles!='')&(nchar(vec_alleles)>max_no_strands)){
#                 for (strand in (max_no_strands+1):nchar(vec_alleles)){
#                     TABLE_INITIAL_CN[paste('Clone_',I_clone,'_strand_',strand,sep='')]  <- 'NA'
#                 }
#             }
#             vec_cols                                        <- which(grepl(paste('Clone_',I_clone,sep=''),colnames(TABLE_INITIAL_CN)))
#             TABLE_INITIAL_CN[vec_rows,vec_cols]             <- 'NA'
#             if (nchar(vec_alleles)==0){
#                 next
#             }
#
#             for (strand in 1:nchar(vec_alleles)){
#
# print(strand)
#
#                 allele                                      <- substr(vec_alleles,strand,strand)
#                 col                                         <- which(colnames(TABLE_INITIAL_CN)==paste('Clone_',I_clone,'_strand_',strand,sep=''))
#                 TABLE_INITIAL_CN[vec_rows,col]              <- allele
#             }
#         }
#     }









#----------------------------------------Output the model variable files
    MODEL_VARIABLES$initial_others                          <- TABLE_INITIAL_OTHERS
    MODEL_VARIABLES$initial_cn                              <- TABLE_INITIAL_CN
    return(MODEL_VARIABLES)
}
