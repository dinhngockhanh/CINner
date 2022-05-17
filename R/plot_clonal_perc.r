plot_clonal_perc <- function(model = "",
                             n_simulations = 0,
                             width = 1000,
                             height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #----------------------------------------Extract true clustering
        phylogeny_clustering_truth <- simulation$sample_phylogeny$phylogeny_clustering_truth
        sample_clustering <- phylogeny_clustering_truth$clustering
        #-------------------------Get list of sample ID's and clone ID's
        sample_clustering$sample_id <- sub("-.*", "", sample_clustering$cell_id)
        list_sample_id <- unique(sample_clustering$sample_id)
        list_clone_id <- unique(sample_clustering$clone_id)
        #------------------------------Find clonal percentages over time
        all_index <- c()
        all_sample_id <- c()
        all_clone_perc <- c()
        all_clone_id <- c()
        for (sample in 1:length(list_sample_id)) {
            sample_id <- list_sample_id[sample]
            clone_perc <- rep(0, length = length(list_clone_id))
            for (clone in 1:length(list_clone_id)) {
                clone_id <- list_clone_id[clone]
                clone_perc[clone] <- length(
                    which(sample_clustering$clone_id == clone_id &
                        sample_clustering$sample_id == sample_id)
                )
            }
            clone_perc <- 100 * clone_perc / sum(clone_perc)

            all_index <- c(all_index, rep(sample, length = length(clone_perc)))
            all_sample_id <- c(all_sample_id, rep(sample_id, length = length(clone_perc)))
            all_clone_perc <- c(all_clone_perc, clone_perc)
            all_clone_id <- c(all_clone_id, list_clone_id)
        }
        table_clone_perc <- data.frame(
            all_index,
            all_sample_id,
            all_clone_perc,
            all_clone_id
        )
        #--------------------------Plot the clonal percentages over time
        filename <- paste(model, "_sim", i, "_clonal_percentages", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)

        print(table_clone_perc)

        p <- ggplot(table_clone_perc, aes(x = all_index, y = all_clone_perc, fill = factor(all_clone_id, levels = list_clone_id))) +
            geom_area(size = 0.5, alpha = 0.8, color = "black") +
            xlab("Samples") +
            ylab("Percentages") +
            labs(fill = "Clone") +
            theme(text = element_text(size = 20)) +
            scale_x_continuous(
                breaks = (1:length(list_sample_id)), labels = list_sample_id, expand = c(0, 0),
                limits = c(1, length(list_sample_id))
            ) +
            scale_y_continuous(
                expand = c(0, 0),
                limits = c(0, 100)
            )
        print(p)
        dev.off()
    }
}
