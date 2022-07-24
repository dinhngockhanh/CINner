#' @export
statistics_hmm <- function(model = "",
                           n_simulations = 0) {
    vec_hmm_ploidy_small <- rep(0, length = n_simulations)
    vec_hmm_ploidy_right <- rep(0, length = n_simulations)
    vec_hmm_ploidy_big <- rep(0, length = n_simulations)
    vec_hmm_cn_small <- rep(0, length = n_simulations)
    vec_hmm_cn_right <- rep(0, length = n_simulations)
    vec_hmm_cn_big <- rep(0, length = n_simulations)
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        cn_profiles_long_TRUTH <- simulation$sample$cn_profiles_long
        cn_profiles_long_HMM <- simulation$sample$cn_profiles_long_hmm
        #---------------------Find ploidy for each cell from TRUTH & HMM
        list_cell_id <- unique(cn_profiles_long_TRUTH$cell_id)
        list_cell_ploidy_TRUTH <- rep(0, length = length(list_cell_id))
        list_cell_ploidy_HMM <- rep(0, length = length(list_cell_id))
        for (j in 1:length(list_cell_id)) {
            cell_id <- list_cell_id[j]
            list_cell_ploidy_TRUTH[j] <- round(mean(cn_profiles_long_TRUTH$state[which(cn_profiles_long_TRUTH$cell_id == cell_id)]))
            list_cell_ploidy_HMM[j] <- round(mean(cn_profiles_long_HMM$state[which(cn_profiles_long_HMM$cell_id == cell_id)]))
        }
        #--------------Find cell counts of correct/wrong inferred ploidy
        hmm_ploidy_small <- length(which((list_cell_ploidy_HMM - list_cell_ploidy_TRUTH) < 0))
        hmm_ploidy_right <- length(which((list_cell_ploidy_HMM - list_cell_ploidy_TRUTH) == 0))
        hmm_ploidy_big <- length(which((list_cell_ploidy_HMM - list_cell_ploidy_TRUTH) > 0))
        #---------------------------Find bin counts of wrong inferred CN
        #   Limit to cells with correct inferred ploidy
        list_cell_id_correct_ploidy <- list_cell_id[which((list_cell_ploidy_HMM - list_cell_ploidy_TRUTH) == 0)]
        mini_TRUTH <- cn_profiles_long_TRUTH[which(cn_profiles_long_TRUTH$cell_id %in% list_cell_id_correct_ploidy), ]
        mini_TRUTH <- mini_TRUTH[order(mini_TRUTH$cell_id, mini_TRUTH$chr, mini_TRUTH$start), ]
        mini_HMM <- cn_profiles_long_HMM[which(cn_profiles_long_HMM$cell_id %in% list_cell_id_correct_ploidy), ]
        mini_HMM <- mini_HMM[order(mini_HMM$cell_id, mini_HMM$chr, mini_HMM$start), ]
        #   Find bin counts of wrong inferred CN per cell with correct inferred ploidy
        hmm_cn_small <- length(which((mini_HMM$state - mini_TRUTH$state) < 0))
        hmm_cn_right <- length(which((mini_HMM$state - mini_TRUTH$state) == 0))
        hmm_cn_big <- length(which((mini_HMM$state - mini_TRUTH$state) > 0))
        #---------------------------Save statistics to simulation record
        simulation$sample$hmm_ploidy_small <- hmm_ploidy_small
        simulation$sample$hmm_ploidy_right <- hmm_ploidy_right
        simulation$sample$hmm_ploidy_big <- hmm_ploidy_big
        simulation$sample$hmm_cn_small <- hmm_cn_small
        simulation$sample$hmm_cn_right <- hmm_cn_right
        simulation$sample$hmm_cn_big <- hmm_cn_big

        vec_hmm_ploidy_small[i] <- hmm_ploidy_small
        vec_hmm_ploidy_right[i] <- hmm_ploidy_right
        vec_hmm_ploidy_big[i] <- hmm_ploidy_big
        vec_hmm_cn_small[i] <- hmm_cn_small
        vec_hmm_cn_right[i] <- hmm_cn_right
        vec_hmm_cn_big[i] <- hmm_cn_big
        #-----------------------------------------Save simulation record
        save(simulation, file = filename)



        print("==============")
        print(i)
        print("--------------")
        print(hmm_ploidy_small)
        print(hmm_ploidy_right)
        print(hmm_ploidy_big)
        print("--------------")
        print(hmm_cn_small)
        print(hmm_cn_right)
        print(hmm_cn_big)
    }
    hmm_stats <- list()

    hmm_stats$ploidy_small <- vec_hmm_ploidy_small
    hmm_stats$ploidy_right <- vec_hmm_ploidy_right
    hmm_stats$ploidy_big <- vec_hmm_ploidy_big
    hmm_stats$cn_small <- vec_hmm_cn_small
    hmm_stats$cn_right <- vec_hmm_cn_right
    hmm_stats$cn_big <- vec_hmm_cn_big

    return(hmm_stats)
}
