#=============================SET ALL PARAMETERS REQUIRED FOR SIMULATION
SIMULATOR_VARIABLES_for_simulation <- function(model) {
#--------------------------------------Set up variables for cancer model
#---Input table of variables from file
    filename                                <- paste(model,'-input-variables.csv',sep='')
    TABLE_VARIABLES                         <- read.table(filename,header=TRUE,sep=',')
#---Input table of total population dynamics
    filename                                <- paste(model,'-input-population-dynamics.csv',sep='')
    TABLE_POPULATION_DYNAMICS               <- read.table(filename,header=TRUE,sep=',')
#---Input table of chromosome bin counts and centromere locations
    filename                                <- paste(model,'-input-copy-number-blocks.csv',sep='')
    TABLE_CHROMOSOME_CN_INFO                <- read.table(filename,header=TRUE,sep=',')
#---Input table of mutational and CNA genes, Warning : Execute later because it is empty now
    filename                                <- paste(model,'-input-cancer-genes.csv',sep='')
    if (file.size(filename)!=0){
        TABLE_CANCER_GENES                  <- read.table(filename,header=TRUE,sep=',')
    }
    else{
        TABLE_CANCER_GENES                  <- data.frame()
    }
#---Input table of CN profiles for the initial population
    filename                                <- paste(model,'-input-initial-cn-profiles.csv',sep='')
    TABLE_INITIAL_COPY_NUMBER_PROFILES      <- read.table(filename,header=TRUE,sep=',')
#---Input other information for the initial population
    filename                                <- paste(model,'-input-initial-others.csv',sep='')
    TABLE_INITIAL_OTHERS                    <- read.table(filename,header=TRUE,sep=',')
#---Set up individual variables from table
    for (i in 1:nrow(TABLE_VARIABLES)) {
        assign(TABLE_VARIABLES[i,1],TABLE_VARIABLES[i,2],envir=.GlobalEnv)
    }
    T_start_time                            <<- age_birth*365
    T_end_time                              <<- age_end*365







#   level_purity is only for SIMULATOR_CLONAL and SIMULATOR_ODE
    level_purity                            <<- 1







#---Set up individual variables from table
    N_chromosomes                           <<- nrow(TABLE_CHROMOSOME_CN_INFO)
    for (i in 1:ncol(TABLE_CHROMOSOME_CN_INFO)) {
        assign(colnames(TABLE_CHROMOSOME_CN_INFO)[i],TABLE_CHROMOSOME_CN_INFO[,i],envir=.GlobalEnv)
    }
    vec_CN_block_no                         <<- Bin_count
    vec_centromere_location                 <<- Centromere_location
#---Set up mutational and CNA driver library (without selection rates)
    if(length(TABLE_CANCER_GENES)==0){
        driver_library                      <<- data.frame()
    }
    else{
        driver_library                      <<- TABLE_CANCER_GENES
    }
#-----------------------------------Set up initial state for simulations
#   Get number of clones in the initial population
    initial_N_clones                        <<- nrow(TABLE_INITIAL_OTHERS)
    vec_header                              <- names(TABLE_INITIAL_COPY_NUMBER_PROFILES)
#---Initialize ID and population for each clone
    initial_clonal_ID                       <<- 1:initial_N_clones
    initial_population                      <<- rep(0,initial_N_clones)
    for (clone in 1:initial_N_clones){
#       Get driver profile of this clone
        loc                                 <- which(TABLE_INITIAL_OTHERS$Clone==clone)
        population                          <- TABLE_INITIAL_OTHERS$Cell_count[loc]
        initial_population[loc]             <- population
    }

#---Initialize the genotypes for each clone
    initial_ploidy_chrom                    <- vector('list',length=initial_N_clones)
    initial_ploidy_allele                   <- vector('list',length=initial_N_clones)
    initial_ploidy_block                    <- vector('list',length=initial_N_clones)
    initial_driver_count                    <- rep(0,initial_N_clones)
    initial_driver_map                      <- vector('list',length=initial_N_clones)
    initial_DNA_length                      <- vector('list',length=initial_N_clones)
    initial_selection_rate                  <- rep(0,initial_N_clones)
    initial_prob_new_drivers                <- rep(0,initial_N_clones)
#---Set up the initial clones' CN genotypes
    for (clone in 1:initial_N_clones){
#       Extract mini table for the CN genotypes of this clone
        text_clone_ID                       <- paste('Clone_',clone,sep='')
        vec_loc                             <- c(1,2,grep(text_clone_ID,vec_header))
        CLONE_INITIAL_COPY_NUMBER_PROFILES  <- TABLE_INITIAL_COPY_NUMBER_PROFILES[,vec_loc]
#       Set up clone's CN genotype
        ploidy_chrom                        <- rep(0,N_chromosomes)
        ploidy_block                        <- list()
        ploidy_allele                       <- list()
        # for (chrom in 1:N_chromosomes){
        for (chrom in 1:1){
#           Get CN genotype for this chromosome
            CHROM_COPY_NUMBER_PROFILES      <- CLONE_INITIAL_COPY_NUMBER_PROFILES[CLONE_INITIAL_COPY_NUMBER_PROFILES$Chromosome==chrom,]
#           Clean CN genotype of unnecessary strands
            vec_delete                      <- c()
            for (column in 3:ncol(CHROM_COPY_NUMBER_PROFILES)){

print(unique(CHROM_COPY_NUMBER_PROFILES[,column]))

                if (all(grepl('NA',CHROM_COPY_NUMBER_PROFILES[,column]))){
                # if (all(is.nan(CHROM_COPY_NUMBER_PROFILES[,c(column)]))){
                    vec_delete              <- c(vec_delete,column)
                }
            }




        }


print(CHROM_COPY_NUMBER_PROFILES)
print(vec_delete)

    }


















#---Set up total population dynamics as function of age (in days)
    for (i in 1:ncol(TABLE_POPULATION_DYNAMICS)) {
        assign(colnames(TABLE_POPULATION_DYNAMICS)[i],TABLE_POPULATION_DYNAMICS[,i],envir=.GlobalEnv)
    }
    vec_age_in_days                         <- 365*Age_in_year
    vec_total_cell_count                    <- Total_cell_count
    linear_app_fun                          <<- approxfun(c(-100*365, vec_age_in_days, 200*365),c(vec_total_cell_count[1], vec_total_cell_count, tail(vec_total_cell_count,1)),method="linear")
    func_expected_population                <<- function(time) linear_app_fun(time)
#---Set up event rate as function of time (in days)
    func_event_rate                         <<- function(time) 1/cell_lifespan
}
