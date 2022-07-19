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
            xlab("GC content") +
            ylab("Readcounts") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 20))
        print(p)
        dev.off()
    }
}
