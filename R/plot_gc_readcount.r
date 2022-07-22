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
<<<<<<< HEAD
<<<<<<< HEAD
        df_plot <- noisy_cn_profiles_long[, c("gc", "map", "reads")]
        vec_delete <- which(df_plot$gc < 0 | df_plot$map < 0)
        if (length(vec_delete) > 0) {
            df_plot <- df_plot[-vec_delete, ]
        }
        #------------------------------Plot readcounts versus GC content
        filename <- paste(model, "_sim", i, "_reads_vs_GC", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- ggplot(df_plot, aes(x = gc, y = reads)) +
            geom_point(colour = "blue", size = 2) +
=======
=======
>>>>>>> 4990b5b (New driver list for HGSOC)



        #-------------------------------Compute cell ploidy = average CN
        noisy_cn_profiles_long$true_ploidy <- 0
        list_cell_id <- unique(noisy_cn_profiles_long$cell_id)
        for (j in 1:length(list_cell_id)) {
            vec_loc <- which(noisy_cn_profiles_long$cell_id == list_cell_id[j])
            noisy_cn_profiles_long$true_ploidy[vec_loc] <- round(mean(noisy_cn_profiles_long$true_CN[vec_loc]))
        }
        noisy_cn_profiles_long$multiplier <- noisy_cn_profiles_long$true_CN / noisy_cn_profiles_long$true_ploidy
        noisy_cn_profiles_long$CN <- paste(noisy_cn_profiles_long$multiplier, "x ploidy", sep = "")



        #--------------------------------Extract data frame for plotting
        df_plot <- noisy_cn_profiles_long[, c("gc", "map", "reads", "CN")]
        df_plot$CN <- factor(df_plot$CN, levels = paste(sort(unique(noisy_cn_profiles_long$multiplier), decreasing = TRUE), "x ploidy", sep = ""))
        vec_delete <- which(df_plot$gc < 0 | df_plot$map < 0 | df_plot$reads == NA | df_plot$reads == NaN | df_plot$reads == Inf)
        if (length(vec_delete) > 0) {
            df_plot <- df_plot[-vec_delete, ]
        }


        #------------------------------Plot readcounts versus GC content
        filename <- paste(model, "_sim", i, "_reads_vs_GC", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- ggplot(df_plot, aes(x = gc, y = reads, col = CN)) +
            geom_point(alpha = 0.05) +
            geom_smooth(na.rm = TRUE, method = "lm", se = TRUE) +
<<<<<<< HEAD
>>>>>>> 4990b5b... New driver list for HGSOC
=======
>>>>>>> 4990b5b (New driver list for HGSOC)
            xlab("GC content") +
            ylab("Readcounts") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20))
        print(p)
        dev.off()
    }
}
