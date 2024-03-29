#' @export
BUILD_initial_population <- function(model_variables = list(),
                                     cell_count = 1,
                                     CN_matrix = list(),
                                     drivers = list()) {
    #------------------------------Update initial clones - other information
    DRIVERS <- ""
    if (length(drivers) > 0) {
        for (i in 1:length(drivers)) {
            strand <- drivers[[i]][[1]]
            unit <- drivers[[i]][[2]]
            driver_ID <- drivers[[i]][[3]]
            if (DRIVERS != "") {
                DRIVERS <- paste(DRIVERS, ";", sep = "")
            }
            DRIVERS <- paste(DRIVERS, driver_ID, "_strand", strand, "_unit", unit, sep = "")
        }
    }
    if (is.null(model_variables$initial_others)) {
        columns <- c("Clone", "Cell_count", "Drivers")

        TABLE_INITIAL_OTHERS <- data.frame(1, cell_count, DRIVERS)
        colnames(TABLE_INITIAL_OTHERS) <- columns
        I_clone <- 1
    } else {
        TABLE_INITIAL_OTHERS <- model_variables$initial_others
        I_clone <- nrow(TABLE_INITIAL_OTHERS) + 1

        TABLE_INITIAL_OTHERS[I_clone, ] <- c(I_clone, cell_count, DRIVERS)
    }
    #------------------------------------Update initial clones - CN profiles
    CN_matrix$Clone <- I_clone
    if (is.null(model_variables$initial_cn)) {
        TABLE_INITIAL_CN <- CN_matrix
    } else {
        TABLE_INITIAL_CN <- model_variables$initial_cn
        TABLE_INITIAL_CN <- rbind(TABLE_INITIAL_CN, CN_matrix)
    }
    #----------------------------------------Output the model variable files
    model_variables$initial_others <- TABLE_INITIAL_OTHERS
    model_variables$initial_cn <- TABLE_INITIAL_CN
    return(model_variables)
}
