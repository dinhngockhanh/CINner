#' @export
plot_gc_readcount <- function(model = "",
                              n_simulations = 0,
                              width = 1000,
                              height = 1000) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #----------------------Input the data frame of GC and readcounts
        noisy_cn_profiles_long <- simulation$sample$noisy_cn_profiles_long
        #-------------------------------Compute cell ploidy = average CN
        noisy_cn_profiles_long$true_ploidy <- 0
        list_cell_id <- unique(noisy_cn_profiles_long$cell_id)
        for (j in 1:length(list_cell_id)) {
            vec_loc <- which(noisy_cn_profiles_long$cell_id == list_cell_id[j])
            noisy_cn_profiles_long$true_ploidy[vec_loc] <- round(mean(noisy_cn_profiles_long$true_CN[vec_loc]))
        }
        noisy_cn_profiles_long$multiplier <- noisy_cn_profiles_long$true_CN / noisy_cn_profiles_long$true_ploidy
        #--------------------------------Extract data frame for plotting
        df_plot <- noisy_cn_profiles_long[, c("gc", "map", "reads", "multiplier")]
        vec_delete <- which(df_plot$gc < 0 | df_plot$map < 0 | df_plot$reads == NA | df_plot$reads == NaN | df_plot$reads == Inf)
        if (length(vec_delete) > 0) {
            df_plot <- df_plot[-vec_delete, ]
        }
        #------------------------------Plot readcounts versus GC content
        df_plot <- df_plot[which(((df_plot$multiplier %% 0.5) == 0) & (df_plot$multiplier <= 1.5)), ]
        df_plot$CN <- factor(paste(df_plot$multiplier, "x ploidy", sep = ""), levels = paste(sort(unique(df_plot$multiplier), decreasing = TRUE), "x ploidy", sep = ""))
        filename <- paste(model, "_sim", i, "_reads_vs_GC", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- ggplot(df_plot, aes(x = gc, y = reads, col = CN)) +
            geom_point(alpha = 0.05) +
            geom_smooth(na.rm = TRUE, method = "lm", se = TRUE) +
            xlab("GC content") +
            ylab("Readcounts") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, 3000), expand = c(0, 0))
        print(p)
        dev.off()
    }
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    # filename <- paste(model_name, "_inferred_ploidy_in_HMM_var=sigma.jpeg", sep = "")
    # jpeg(file = filename, width = 2000, height = 1000)
    df<-data.frame(
        Sigma = rep(vec.sigma1[1:5], 2),
        Percentage = c(
            100 * hmm_stats$ploidy_small[1:5] / (hmm_stats$ploidy_small[1:5] + hmm_stats$ploidy_right[1:5] + hmm_stats$ploidy_big[1:5]),
            100 * hmm_stats$ploidy_big[1:5] / (hmm_stats$ploidy_small[1:5] + hmm_stats$ploidy_right[1:5] + hmm_stats$ploidy_big[1:5])
        ),
        Group = c(rep("Inferred ploidy is too small", 5), rep("Inferred ploidy is too big", 5))
    )
    df$Sigma <- factor(df$Sigma, levels=vec.sigma1[1:5])
    p <- ggplot(
        data = df,
        aes(x = Sigma, y = Percentage, fill = Group)
    ) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20), legend.position = "bottom")

    # print(p)
    # dev.off()
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
    ############################################################
}
