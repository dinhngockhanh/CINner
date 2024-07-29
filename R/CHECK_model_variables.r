#' Check/Clean up the model variables
#' 
#' @description
#' `CHECK_model_variable` returns an updated object that checks the specified model variables to clean up the libraries of CN information and chromosome arm selection rates.
#'
#' @inheritParams BUILD_initial_population
#' 
#' @examples 
#' 
#' CHECK_model_variables(model_variables)
#'
#'
#' @export
CHECK_model_variables <- function(model_variables) {
    TABLE_INITIAL_CN <- model_variables$initial_cn
    #---------------------------------Clean up library of CN information
    TABLE_CHROMOSOME_CN_INFO <- model_variables$cn_info
    # print(TABLE_INITIAL_CN)
    # print(TABLE_CHROMOSOME_CN_INFO)
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
    return(model_variables)
}
