#' @export
SIMULATOR_FULL_PHASE_1_clonal_population_cleaning <- function() {
    vec_delete <- which(clonal_population_next == 0)
    if (length(vec_delete) != 0) {
        clonal_population_current <<- clonal_population_current[-vec_delete]
        clonal_population_next <<- clonal_population_next[-vec_delete]
        clonal_ID_current <<- clonal_ID_current[-vec_delete]
    }
}
