SAVE_model_variables <- function(model_name='',
                                model_variables=list()){


    TABLE_CHROMOSOME_CN_INFO        <- model_variables$cn_info
    TABLE_INITIAL_CN                <- model_variables$initial_cn

    vec_delete                      <- c()
    for (i in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                       <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i]
        vec_loc                     <- which(TABLE_INITIAL_CN$Chromosome==chrom)
        if (length(vec_loc)==0){
            vec_delete              <- c(vec_delete,chrom)
        }
    }

    if (length(vec_delete)>0){
        model_variables$cn_info     <- model_variables$cn_info[-which(is.element(model_variables$cn_info$Chrom,vec_delete)),]
    }

#---Save file for general variables
    filename    <- paste(model_name,'-input-variables.csv',sep='')
    write.csv(model_variables$general_variables,filename,row.names=FALSE)
#---Save file for CN information
    filename    <- paste(model_name,'-input-copy-number-blocks.csv',sep='')
    write.csv(model_variables$cn_info,filename,row.names=FALSE)
#---Save file for population dynamics
    filename    <- paste(model_name,'-input-population-dynamics.csv',sep='')
    write.csv(model_variables$population_dynamics,filename,row.names=FALSE)
#---Save file for sampling information
    filename    <- paste(model_name,'-input-sampling.csv',sep='')
    write.csv(model_variables$sampling_info,filename,row.names=FALSE)
#---Save file for driver library
    filename    <- paste(model_name,'-input-cancer-genes.csv',sep='')
    write.csv(model_variables$driver_library,filename,row.names=FALSE)
#---Save file for initial clones' CN profiles
    filename    <- paste(model_name,'-input-initial-cn-profiles.csv',sep='')
    write.csv(model_variables$initial_cn,filename,row.names=FALSE)
#---Save file for initial clones' other information
    filename    <- paste(model_name,'-input-initial-others.csv',sep='')
    write.csv(model_variables$initial_others,filename,row.names=FALSE)
}
