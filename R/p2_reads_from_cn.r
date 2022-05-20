p2_reads_from_cn <- function(simulation) {
    noisy_cn_profiles_long <- simulation$sample$sample_genotype_profiles
    #-------------------Find true_CN and BAF for every bin in every cell
    #   Find true_CN and BAF for every bin in every cell
    noisy_cn_profiles_long$true_CN <- noisy_cn_profiles_long$state
    noisy_cn_profiles_long$BAF <- noisy_cn_profiles_long$Min / noisy_cn_profiles_long$state
    #   Eliminate unnecessary columns
    noisy_cn_profiles_long <- noisy_cn_profiles_long[c("chr", "start", "end", "cell_id", "true_CN", "BAF")]
    #-----------------------------------------Find list of all cell ID's
    all_cell_id <- unique(noisy_cn_profiles_long$cell_id)
    #--------------------Find GC/mappability for every bin in every cell
    cat("A===Simulate noisy CN with GC/mappability biases...\n")
    cat("+---Assign GC to every bin in every cell...\n")
    #   Find GC/mappability for the first cell
    cell_id <- all_cell_id[1]
    one_noisy_cn_profiles_long <- data.frame(
        chr = noisy_cn_profiles_long$chr[noisy_cn_profiles_long$cell_id == cell_id],
        start = noisy_cn_profiles_long$start[noisy_cn_profiles_long$cell_id == cell_id],
        end = noisy_cn_profiles_long$end[noisy_cn_profiles_long$cell_id == cell_id]
    )
    one_noisy_cn_profiles_long$gc <- -1
    one_noisy_cn_profiles_long$map <- -1
    for (row in 1:nrow(one_noisy_cn_profiles_long)) {
        chr <- one_noisy_cn_profiles_long$chr[row]
        start <- one_noisy_cn_profiles_long$start[row]
        end <- one_noisy_cn_profiles_long$end[row]
        loc <- which(gc$chr == chr & gc$start == start & gc$end == end)
        if (length(loc) > 0) {
            one_noisy_cn_profiles_long$gc[row] <- gc$gc[loc]
            one_noisy_cn_profiles_long$map[row] <- gc$map[loc]
        }
    }
    #   Assign GC/mappability for every cell
    noisy_cn_profiles_long$gc <- -1
    noisy_cn_profiles_long$map <- -1
    for (cell in 1:length(all_cell_id)) {
        cell_id <- all_cell_id[cell]
        vec_loc <- which(noisy_cn_profiles_long$cell_id == cell_id)
        noisy_cn_profiles_long$gc[vec_loc] <- one_noisy_cn_profiles_long$gc
        noisy_cn_profiles_long$map[vec_loc] <- one_noisy_cn_profiles_long$map
    }
    #   Delete CN bins without GC/mappability information
    vec_delete <- which(noisy_cn_profiles_long$gc < 0)
    noisy_cn_profiles_long <- noisy_cn_profiles_long[-vec_delete, ]
    # #--------------------Find GC/mappability for every bin in every cell
    # #   Find GC/mappability for every bin in every cell
    # noisy_cn_profiles_long$gc <- -1
    # noisy_cn_profiles_long$map <- -1
    # cat("+---Assign GC to every bin in every cell...\n")
    # pb <- txtProgressBar(
    #     min = 0, max = nrow(gc),
    #     style = 3, width = 50, char = "="
    # )
    # for (row in 1:nrow(gc)) {
    #     setTxtProgressBar(pb, row)
    #
    #     vec_loc <- which(noisy_cn_profiles_long$chr == gc$chr[row] &
    #         noisy_cn_profiles_long$start == gc$start[row] &
    #         noisy_cn_profiles_long$end == gc$end[row])
    #     if (length(vec_loc) > 0) {
    #         noisy_cn_profiles_long$gc[vec_loc] <- gc$gc[row]
    #         noisy_cn_profiles_long$map[vec_loc] <- gc$map[row]
    #     }
    # }
    # cat("\n")
    # #   Delete CN bins without GC/mappability information
    # vec_delete <- which(noisy_cn_profiles_long$gc < 0)
    # noisy_cn_profiles_long <- noisy_cn_profiles_long[-vec_delete, ]
    #-------------Simulate noisy CN profiles for each cell in the sample
    #-------------------------------------with GC and mappability biases
    #   Model GC bias
    cat("+---Model GC bias...\n")
    noisy_cn_profiles_long$observed_CN <- noisy_cn_profiles_long$true_CN *
        (gc_slope * noisy_cn_profiles_long$gc + gc_int)
    #   Model random noise in observed CN
    cat("+---Model random noise in observed CN...\n")
    noisy_cn_profiles_long$noisy_CN <- rgamma(
        n = nrow(noisy_cn_profiles_long),
        shape = noisy_cn_profiles_long$observed_CN / sigma1,
        scale = sigma1
    )
    #   Normalize noisy CN per cell
    cat("+---Normalize noisy CN per cell...\n")
    noisy_cn_profiles_long$noisy_CN_pval <- 0
    for (row in 1:length(all_cell_id)) {
        vec_loc <- which(noisy_cn_profiles_long$cell_id == all_cell_id[row])
        noisy_cn_profiles_long$noisy_CN_pval[vec_loc] <-
            noisy_cn_profiles_long$noisy_CN[vec_loc] / sum(noisy_cn_profiles_long$noisy_CN[vec_loc])
    }
    #   Model total read counts per bin per cell
    cat("+---Model total read counts...\n")
    noisy_cn_profiles_long$reads <- 0
    for (row in 1:length(all_cell_id)) {
        vec_loc <- which(noisy_cn_profiles_long$cell_id == all_cell_id[row])
        noisy_cn_profiles_long$reads[vec_loc] <- rmultinom(
            n = 1,
            size = num_reads,
            prob = noisy_cn_profiles_long$noisy_CN_pval[vec_loc]
        )
    }
    #   Model minor read counts per bin per cell
    cat("+---Model minor read counts...\n")
    noisy_cn_profiles_long$minor_reads <- rbinom(
        n = nrow(noisy_cn_profiles_long),
        size = noisy_cn_profiles_long$reads,
        prob = noisy_cn_profiles_long$BAF
    )
    #   Find major read counts per bin per cell
    cat("+---Find major read counts...\n")
    noisy_cn_profiles_long$major_reads <-
        noisy_cn_profiles_long$reads -
        noisy_cn_profiles_long$minor_reads
    #----------------------------------Use HMMcopy to find inferred copy
    cat("B===Use HMMcopy to infer bias-free CN...\n")
    noisy_cn_profiles_long_list <- vector("list", length = length(all_cell_id))
    for (cell in 1:length(all_cell_id)) {
        #---Get noisy CN profile for this cell
        cell_id <- all_cell_id[cell]
        vec_loc <- which(noisy_cn_profiles_long$cell_id == cell_id)
        uncorrected_reads <- noisy_cn_profiles_long[vec_loc, ]
        uncorrected_reads <- uncorrected_reads[c("chr", "start", "end", "reads", "gc", "map", "true_CN")]
        #---Correct noisy CN profile for GC and mappability biases
        corrected_copy <- correctReadcount(uncorrected_reads)
        #---Run HMMcopy to infer CN state
        #   Prepare parameters for (modified) HMMcopy pipeline
        cell <- cell_id

        param <- HMMsegment(corrected_copy, getparam = TRUE)
        class(param) <- "data.table"

        multipliers <- "2"


        # opt <- list()
        # opt$param_str <- 2
        # # ???????????????????????????
        # opt$param_e <- 0.9
        # # opt$param_e <- 2
        # opt$param_mu <- "2"
        # opt$param_l <- 2
        # opt$param_nu <- 2
        # # ???????????????????????????
        # opt$param_k <- "1"
        # # opt$param_k <- "2"
        # opt$param_m <- "2"
        # opt$param_eta <- "2"
        # opt$param_g <- 2
        # opt$param_s <- 2
        # param <- get_parameters_MODIFIED(opt$param_str, opt$param_e, opt$param_mu, opt$param_l, opt$param_nu, opt$param_k, opt$param_m, opt$param_eta, opt$param_g, opt$param_s)


        #   Run (modified) HMMcopy pipeline
        output <- run_hmmcopy_MODIFIED(cell, corrected_copy, param, multipliers)

        # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # str(output)
        # print(output$auto_ploidy.reads)

        #   Put the data together for this cell
        corrected_copy <- corrected_copy[order(corrected_copy$chr, corrected_copy$start), ]
        colnames(corrected_copy)[which(names(corrected_copy) == "copy")] <- "copy_log2"
        corrected_copy$multiplier <- output$auto_ploidy.reads$multiplier
        corrected_copy$copy <- output$auto_ploidy.reads$copy
        corrected_copy$state <- output$auto_ploidy.reads$state
        #   Store inferred CN profiles for this cell
        noisy_cn_profiles_long_list[[cell]] <- corrected_copy
    }
    noisy_cn_profiles_long <- rbindlist(noisy_cn_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
    #-------------------------------------------Output noisy CN profiles
    simulation$sample$noisy_cn_profiles_long <- noisy_cn_profiles_long
    return(simulation)
}
