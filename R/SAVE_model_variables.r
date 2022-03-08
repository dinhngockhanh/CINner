SAVE_model_variables <- function(MODEL_NAME='',
                                MODEL_VARIABLES=list()){


    TABLE_CHROMOSOME_CN_INFO        <- MODEL_VARIABLES$cn_info
    TABLE_INITIAL_CN                <- MODEL_VARIABLES$initial_cn

    vec_delete                      <- c()
    for (i in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                       <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i]
        vec_loc                     <- which(TABLE_INITIAL_CN$Chromosome==chrom)
        if (length(vec_loc)==0){
            vec_delete              <- c(vec_delete,chrom)
        }
    }

    if (length(vec_delete)>0){
        MODEL_VARIABLES$cn_info     <- MODEL_VARIABLES$cn_info[-which(is.element(MODEL_VARIABLES$cn_info$Chrom,vec_delete)),]
    }

#---Save file for general variables
    filename    <- paste(MODEL_NAME,'-input-variables.csv',sep='')
    write.csv(MODEL_VARIABLES$general_variables,filename,row.names=FALSE)
#---Save file for CN information
    filename    <- paste(MODEL_NAME,'-input-copy-number-blocks.csv',sep='')
    write.csv(MODEL_VARIABLES$cn_info,filename,row.names=FALSE)
#---Save file for population dynamics
    filename    <- paste(MODEL_NAME,'-input-population-dynamics.csv',sep='')
    write.csv(MODEL_VARIABLES$population_dynamics,filename,row.names=FALSE)
#---Save file for driver library
    filename    <- paste(MODEL_NAME,'-input-cancer-genes.csv',sep='')
    write.csv(MODEL_VARIABLES$driver_library,filename,row.names=FALSE)
#---Save file for initial clones' CN profiles
    filename    <- paste(MODEL_NAME,'-input-initial-cn-profiles.csv',sep='')
    write.csv(MODEL_VARIABLES$initial_cn,filename,row.names=FALSE)
#---Save file for initial clones' other information
    filename    <- paste(MODEL_NAME,'-input-initial-others.csv',sep='')
    write.csv(MODEL_VARIABLES$initial_others,filename,row.names=FALSE)
}
