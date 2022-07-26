#' @export
SAVE_model_variables <- function(model_name = "",
                                 model_variables = list()) {
    TABLE_INITIAL_CN <- model_variables$initial_cn
    #---------------------------------Clean up library of CN information
    TABLE_CHROMOSOME_CN_INFO <- model_variables$cn_info
    vec_delete <- c()
    for (i in 1:nrow(TABLE_CHROMOSOME_CN_INFO)) {
        chrom <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i]
        vec_loc <- which(TABLE_INITIAL_CN$Chromosome == chrom)
        if (length(vec_loc) == 0) {
            vec_delete <- c(vec_delete, chrom)
        }
    }
    if (length(vec_delete) > 0) {
        model_variables$cn_info <- model_variables$cn_info[-which(is.element(model_variables$cn_info$Chrom, vec_delete)), ]
    }
    #-----------------Clean up library of chromosome arm selection rates
    TABLE_CANCER_ARMS <- model_variables$chromosome_arm_library
    if (!is.null(TABLE_CANCER_ARMS)) {
        vec_delete <- c()
        for (i in 1:nrow(TABLE_CANCER_ARMS)) {
            chrom <- TABLE_CANCER_ARMS$Chromosome[i]
            vec_loc <- which(TABLE_INITIAL_CN$Chromosome == chrom)
            if (length(vec_loc) == 0) {
                vec_delete <- c(vec_delete, chrom)
            }
        }
        if (length(vec_delete) > 0) {
            model_variables$chromosome_arm_library <- model_variables$chromosome_arm_library[-which(is.element(model_variables$chromosome_arm_library$Chromosome, vec_delete)), ]
        }
    }
    #---Save file for general variables
    filename <- paste(model_name, "-input-variables.csv", sep = "")
    write.csv(model_variables$general_variables, filename, row.names = FALSE)
    #---Save file for selection model
    filename <- paste(model_name, "-input-selection-model.csv", sep = "")
    write.csv(model_variables$selection_model, filename, row.names = FALSE)
    #---Save file for GC content and mappability per CN bin
    filename <- paste(model_name, "-input-gc.csv", sep = "")
    write.csv(model_variables$gc_and_mappability, filename, row.names = FALSE)
    #---Save file for CN information
    filename <- paste(model_name, "-input-copy-number-blocks.csv", sep = "")
    write.csv(model_variables$cn_info, filename, row.names = FALSE)
    #---Save file for population dynamics
    filename <- paste(model_name, "-input-population-dynamics.csv", sep = "")
    write.csv(model_variables$population_dynamics, filename, row.names = FALSE)
    #---Save file for sampling information
    filename <- paste(model_name, "-input-sampling.csv", sep = "")
    write.csv(model_variables$sampling_info, filename, row.names = FALSE)
    #---Save file for driver library
    if (!is.null(model_variables$driver_library)) {
        filename <- paste(model_name, "-input-selection-genes.csv", sep = "")
        write.csv(model_variables$driver_library, filename, row.names = FALSE)
    }


    #---Save file for chromosome arm selection library
    if (!is.null(model_variables$chromosome_arm_library)) {
        filename <- paste(model_name, "-input-selection-chrom-arm.csv", sep = "")
        write.csv(model_variables$chromosome_arm_library, filename, row.names = FALSE)
    }



    #---Save file for initial clones' CN profiles
    filename <- paste(model_name, "-input-initial-cn-profiles.csv", sep = "")
    write.csv(model_variables$initial_cn, filename, row.names = FALSE)
    #---Save file for initial clones' other information
    filename <- paste(model_name, "-input-initial-others.csv", sep = "")
    write.csv(model_variables$initial_others, filename, row.names = FALSE)
}
