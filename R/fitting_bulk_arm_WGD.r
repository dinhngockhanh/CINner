#---------------------Function for WGD rate & aneuploidy score from data
get_WGD_stats_from_data <- function(copynumber_DATA,
                                    DATA_wgd,
                                    copynumber_coordinates,
                                    cn_info,
                                    list_targets) {
    CNbin_length <- copynumber_coordinates$end[1] - copynumber_coordinates$start[1] + 1
    #---Find CN profile of every sample (normalized and calibrated for arm level)
    DATA_gainloss <- gainloss_DATA(
        copynumber_DATA, copynumber_coordinates,
        ploidy_normalization = TRUE,
        use_rbindlist = TRUE,
        arm_level = TRUE,
        round = FALSE,
        pos_centromeres = cn_info
    )
    copynumber_data <- DATA_gainloss$copynumber_data
    list_samples <- colnames(copynumber_data)
    list_samples <- list_samples[-c(1:4)]
    df_sample_stats <- data.frame(sample = list_samples)
    #---Find WGD status of every sample
    df_sample_stats$WGD <- 0
    list_delete <- c()
    for (i in 1:nrow(df_sample_stats)) {
        sample <- df_sample_stats$sample[i]
        wgd_uncertain <- DATA_wgd$wgd_uncertain[which(DATA_wgd$samplename == sample)]
        wgd_status <- DATA_wgd$wgd_status[which(DATA_wgd$samplename == sample)]
        if (wgd_uncertain == TRUE) list_delete <- c(list_delete, i)
        if (wgd_status == "wgd") {
            df_sample_stats$WGD[i] <- 1
        } else if (wgd_status == "no_wgd") {
            df_sample_stats$WGD[i] <- 0
        }
    }
    if (length(list_delete) > 0) df_sample_stats <- df_sample_stats[-list_delete, ]
    #---Find FGA of every sample
    df_sample_stats$FGA <- 0
    for (i in 1:nrow(df_sample_stats)) {
        sample <- df_sample_stats$sample[i]
        sample_CN <- copynumber_data[[sample]]
        df_sample_stats$FGA[i] <- length(which(sample_CN != 2)) / length(sample_CN)
    }
    #---Extract statistics
    stat <- list()
    for (i in 1:length(list_targets)) {
        target <- list_targets[i]
        if (target == "WGD_proportion") {
            stat$n_total <- length(df_sample_stats$WGD)
            stat$n_WGD <- length(which(df_sample_stats$WGD == 1))
            stat$n_nonWGD <- length(which(df_sample_stats$WGD == 0))
            stat$WGD_proportion <- length(which(df_sample_stats$WGD == 1)) / length(df_sample_stats$WGD)
        } else if (target == "FGA_difference") {
            if (length(which(df_sample_stats$WGD == 1)) > 0) {
                FGA_WGD <- df_sample_stats$FGA[which(df_sample_stats$WGD == 1)]
            } else {
                FGA_WGD <- 0
            }
            if (length(which(df_sample_stats$WGD == 0)) > 0) {
                FGA_nonWGD <- df_sample_stats$FGA[which(df_sample_stats$WGD == 0)]
            } else {
                FGA_nonWGD <- 0
            }
            stat$FGA_difference <- mean(FGA_WGD) - mean(FGA_nonWGD)
        } else if (target == "WGD_FGA_by_sample") {
            stat$WGD_FGA_by_sample <- df_sample_stats
        }
    }
    #---Return statistics for sample cohort
    return(stat)
}
