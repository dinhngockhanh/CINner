p2_reads_from_cn <- function(simulation, report_progress) {
    noisy_cn_profiles_long <- simulation$sample$cn_profiles_long
    #-------------------Find true_CN and BAF for every bin in every cell
    #   Find true_CN and BAF for every bin in every cell
    noisy_cn_profiles_long$true_CN <- noisy_cn_profiles_long$state
    noisy_cn_profiles_long$true_BAF <- noisy_cn_profiles_long$Min / noisy_cn_profiles_long$state
    #   Eliminate unnecessary columns
    noisy_cn_profiles_long <- noisy_cn_profiles_long[c("chr", "start", "end", "cell_id", "true_CN", "true_BAF")]
    #-----------------------------------------Find list of all cell ID's
    all_cell_id <- unique(noisy_cn_profiles_long$cell_id)
    #--------------------Find GC/mappability for every bin in every cell
    if (report_progress == TRUE) {
        cat("+---Assign GC to every bin in every cell...\n")
    }
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
        loc <- which(table_gc$chr == chr & table_gc$start == start & table_gc$end == end)
        if (length(loc) > 0) {
            one_noisy_cn_profiles_long$gc[row] <- table_gc$gc[loc]
            one_noisy_cn_profiles_long$map[row] <- table_gc$map[loc]
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
    #     vec_loc <- which(noisy_cn_profiles_long$chr == table_gc$chr[row] &
    #         noisy_cn_profiles_long$start == table_gc$start[row] &
    #         noisy_cn_profiles_long$end == table_gc$end[row])
    #     if (length(vec_loc) > 0) {
    #         noisy_cn_profiles_long$gc[vec_loc] <- table_gc$gc[row]
    #         noisy_cn_profiles_long$map[vec_loc] <- table_gc$map[row]
    #     }
    # }
    # cat("\n")
    # #   Delete CN bins without GC/mappability information
    # vec_delete <- which(noisy_cn_profiles_long$gc < 0)
    # noisy_cn_profiles_long <- noisy_cn_profiles_long[-vec_delete, ]
    #----------------------------Find locations with positive GC content
    vec_work <- which(noisy_cn_profiles_long$gc > 0)
    #-------------Simulate noisy CN profiles for each cell in the sample
    #-------------------------------------with GC and mappability biases
    #---Model GC bias
    if (report_progress == TRUE) {
        cat("+---Model GC bias...\n")
    }
    noisy_cn_profiles_long$observed_CN <- -1
    noisy_cn_profiles_long$observed_CN[vec_work] <- noisy_cn_profiles_long$true_CN[vec_work] *
        (gc_slope * noisy_cn_profiles_long$gc[vec_work] + gc_int)
    #---Model random noise in observed CN
    if (report_progress == TRUE) {
        cat("+---Model random noise in observed CN...\n")
    }
    noisy_cn_profiles_long$noisy_CN <- -1
    noisy_cn_profiles_long$noisy_CN[vec_work] <- rgamma(
        n = length(vec_work),
        shape = noisy_cn_profiles_long$observed_CN[vec_work] / sigma1,
        scale = sigma1
    )
    #---Normalize noisy CN per cell
    if (report_progress == TRUE) {
        cat("+---Normalize noisy CN per cell...\n")
    }
    noisy_cn_profiles_long$noisy_CN_pval <- 0
    for (row in 1:length(all_cell_id)) {
        vec_loc <- which(noisy_cn_profiles_long$cell_id == all_cell_id[row] & noisy_cn_profiles_long$gc > 0)
        noisy_cn_profiles_long$noisy_CN_pval[vec_loc] <-
            noisy_cn_profiles_long$noisy_CN[vec_loc] / sum(noisy_cn_profiles_long$noisy_CN[vec_loc])
    }
    #---Model total read counts per bin per cell
    if (report_progress == TRUE) {
        cat("+---Model total read counts...\n")
    }
    noisy_cn_profiles_long$reads <- 0
    for (row in 1:length(all_cell_id)) {
        vec_loc <- which(noisy_cn_profiles_long$cell_id == all_cell_id[row] & noisy_cn_profiles_long$gc > 0)
        noisy_cn_profiles_long$reads[vec_loc] <- rmultinom(
            n = 1,
            size = num_reads,
            prob = noisy_cn_profiles_long$noisy_CN_pval[vec_loc]
        )
    }
    #---Model minor read counts per bin per cell
    if (report_progress == TRUE) {
        cat("+---Model minor read counts...\n")
    }
    noisy_cn_profiles_long$minor_reads <- 0
    noisy_cn_profiles_long$minor_reads[vec_work] <- rbinom(
        n = length(vec_work),
        size = noisy_cn_profiles_long$reads[vec_work],
        prob = noisy_cn_profiles_long$true_BAF[vec_work]
    )
    #---Find major read counts per bin per cell
    if (report_progress == TRUE) {
        cat("+---Find major read counts...\n")
    }
    noisy_cn_profiles_long$major_reads <- 0
    noisy_cn_profiles_long$major_reads[vec_work] <-
        noisy_cn_profiles_long$reads[vec_work] -
        noisy_cn_profiles_long$minor_reads[vec_work]
    #-------------------------------------------Output noisy CN profiles
    simulation$sample$noisy_cn_profiles_long <- noisy_cn_profiles_long
    return(simulation)
}
