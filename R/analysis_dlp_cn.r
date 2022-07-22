analysis_dlp_cn <- function(CNbins, gc) {
    #------------------------------------Get list of cell ID's in sample
    list_cell_id <- unique(CNbins$cell_id)
    #-------------------------Supplement the CN data with GC/mappability
    CNbins_gc <- CNbins
    CNbins_gc$gc <- -1
    CNbins_gc$map <- -1
    for (row in 1:nrow(gc)) {
        vec_loc <- which(CNbins_gc$chr == gc$chr[row] & CNbins_gc$start_point == gc$start[row] & CNbins_gc$end_point == gc$end[row])
        CNbins_gc$gc[vec_loc] <- gc$gc[row]
        CNbins_gc$map[vec_loc] <- gc$map[row]
    }
    vec_delete <- which(CNbins_gc$gc <= 0 | CNbins_gc$map <= 0 | CNbins_gc$state <= 0)
    if (length(vec_delete) > 0) {
        CNbins_gc <- CNbins_gc[-vec_delete, ]
    }
    #-----------------------------------Compute cell ploidy = average CN
    CNbins_gc$true_ploidy <- 0
    for (j in 1:length(list_cell_id)) {
        vec_loc <- which(CNbins_gc$cell_id == list_cell_id[j])
        CNbins_gc$true_ploidy[vec_loc] <- round(mean(CNbins_gc$state[vec_loc]))
    }
    CNbins_gc$multiplier <- CNbins_gc$state / CNbins_gc$true_ploidy

    CNbins_gc <- CNbins_gc[-which(CNbins_gc$multiplier %% 0.5 != 0), ]

    CNbins_gc$CN <- paste(CNbins_gc$multiplier, "x ploidy", sep = "")




    #-------------------
    ggplot(CNbins_gc, aes(x = gc, y = reads, col = CN)) +
        geom_point(alpha = 0.5) +
        geom_smooth(na.rm = TRUE, method = "lm", se = TRUE)









    # for (i in 1:length(list_cell_id)) {
    #     cell_id <- list_cell_id[i]
    #
    #
    #
    #
    #
    #     CNbins_cell <- gc
    #     CNbins_cell$reads <- 0
    #     CNbins_cell$state <- 0
    #     for (row in 1:nrow(CNbins_cell)) {
    #         loc <- which(CNbins$cell_id == cell_id & CNbins$chr == CNbins_cell$chr[row] & CNbins$start_point == CNbins_cell$start[row] & CNbins$end_point == CNbins_cell$end[row])
    #         CNbins_cell$reads[row] <- CNbins$reads[loc]
    #         CNbins_cell$state[row] <- CNbins$state[loc]
    #     }
    #     list_CNbins_cell[[i]] <- CNbins_cell
    #     print(CNbins_cell)
    # }
    # #################
    # #################
    # #################
    # #################
    # #################
    # list_cell_id <<- list_cell_id
    # list_CNbins_cell <<- list_CNbins_cell
    # #################
    # #################
    # #################
    # #################
    # #################
    #
    #
    #
    # # #--------------------------Find distribution of total reads per cell
    # # list_total_reads <- rep(0,length=length(list_cell_id))
    # # for (i in 1:length(list_cell_id)){
    # #     list_total_reads[i]<-sum(CNbins$reads[CNbins$])
    # # }
    # #
    # #
    # #
    # #
    # # print(list_cell_id)
}
