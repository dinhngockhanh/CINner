BUILD_simulator_variables <- function(cell_lifespan     = 4,
                                      T_0               = list(0,'year'),
                                      T_end             = list(100,'year'),
                                      Population_end    = Inf,
                                      Max_events        = Inf,){

print(cell_lifespan)
print(T_0)
print(T_end)

#---------------------------Build model input file for general variables
    columns                         <- c('Variable','Value','Unit','Note')
    TABLE_VARIABLES                 <- data.frame(matrix(nrow=0,ncol=length(columns)))
    colnames(TABLE_VARIABLES)       <- columns
    N_row                           <- 0
#   Set up the start time of simulations
    age_birth                       <- T_0[[1]]
    age_birth_unit                  <- T_0[[2]]
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('age_birth',age_birth,age_birth_unit,'Age when simulation starts')
    if (age_birth_unit=='day'){
        T_start_time                <- age_birth
    }else{if (age_birth_unit=='week'){
        T_start_time                <- 7*age_birth
    }else{if (age_birth_unit=='month'){
        T_start_time                <- 30*age_birth
    }else{if (age_birth_unit=='year'){
        T_start_time                <- 365*age_birth
    }}}}
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('T_start_time',T_start_time,'day','Age when simulation starts (for internal use)')
#   Set up the end time of simulations
    age_end                         <- T_end[[1]]
    age_end_unit                    <- T_end[[2]]
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('age_end',age_end,age_end_unit,'Age when simulation stops')
    if (age_end_unit=='day'){
        T_end_time                  <- age_end
    }else{if (age_end_unit=='week'){
        T_end_time                  <- 7*age_end
    }else{if (age_end_unit=='month'){
        T_end_time                  <- 30*age_end
    }else{if (age_end_unit=='year'){
        T_end_time                  <- 365*age_end
    }}}}
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('T_end_time',T_end_time,'day','Age when simulation stops (for internal use)')
#   Set up condition to end simulation prematurely based on population size or event count
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('Population_end',Population_end,'cell count','Condition for ending simulation (Inf if no condition)')
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('Max_events',Max_events,'event count','Condition for ending simulation (Inf if no condition)')
#   Set up cell lifespan
    N_row                           <- N_row+1
    TABLE_VARIABLES[N_row,]         <- c('cell_lifespan',cell_lifespan,'day','Lifespan of one cell')



    print(TABLE_VARIABLES)
}
