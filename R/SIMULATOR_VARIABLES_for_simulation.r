#' Input variables to simulate a cancer model
#'
#' @param model Name of model.
#' @return Nothing.
#' @examples
#' SIMULATOR_VARIABLES_for_simulation("FALLOPIAN-TUBES-NEUTRAL")
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
    if (is.na(driver_order)){
        driver_order                        <<- list()
    }
    else{
#       ...
#       ...
#       ...
#       ...
#       ...
#       ...
#       ...
#       ...
#       ...
#       ...
    }
    if(length(TABLE_CANCER_GENES)==0){
        driver_library                      <<- list()
    }
    else{
        driver_library                      <<- TABLE_CANCER_GENES
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
