plot_comparison <- function(simulation_statistics_all = list(),
                            vec_time_points = c(0),
                            legend_batches = c(0),
                            plotfile_prefix = "") {
    N_batches <- length(simulation_statistics_all)
    #----------------------------------------Extract values for plotting
    #   Extract clone counts in total population
    df_n_clones_true_in_tot_pop <- data.frame(
        time = c(),
        Parameter = c(),
        true_n_clones_whole_pop_mean = c(),
        true_n_clones_whole_pop_sd = c()
    )
    for (batch in 1:N_batches) {
        batch_statistics <- simulation_statistics_all[[batch]]
        true_n_clones_whole_pop_times <- batch_statistics$true_n_clones_whole_pop_times
        true_n_clones_whole_pop_mean <- batch_statistics$true_n_clones_whole_pop_mean
        true_n_clones_whole_pop_sd <- batch_statistics$true_n_clones_whole_pop_sd
        vec_mean <- rep(0, length = length(vec_time_points))
        vec_sd <- rep(0, length = length(vec_time_points))
        for (pos in 1:length(vec_time_points)) {
            loc <- which.min(abs(true_n_clones_whole_pop_times - vec_time_points[pos]))
            vec_mean[pos] <- true_n_clones_whole_pop_mean[loc]
            vec_sd[pos] <- true_n_clones_whole_pop_sd[loc]
        }
        df_n_clones_true_in_tot_pop_batch <- data.frame(
            time = vec_time_points,
            Parameter = legend_batches[batch],
            true_n_clones_whole_pop_mean = vec_mean,
            true_n_clones_whole_pop_sd = vec_sd
        )
        df_n_clones_true_in_tot_pop <- rbind(df_n_clones_true_in_tot_pop, df_n_clones_true_in_tot_pop_batch)
    }
    #   Extract clone counts in sample
    df_n_clones_in_sample <- data.frame(
        Statistics = c(),
        batch = c(),
        mean = c(),
        sd = c()
    )
    for (batch in 1:N_batches) {
        batch_statistics <- simulation_statistics_all[[batch]]
        df_n_clones_in_sample_batch <- data.frame(
            Statistics = c("Total-CN-based clone count", "True clone count"),
            batch = rep(legend_batches[batch], length = 2),
            mean = c(
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "total_cn_based_n_clones_mean"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "true_n_clones_mean"), 1]
            ),
            sd = c(
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "total_cn_based_n_clones_sd"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "true_n_clones_sd"), 1]
            )
        )
        df_n_clones_in_sample <- rbind(df_n_clones_in_sample, df_n_clones_in_sample_batch)
    }
    #   Extract percentages of each event class in sample
    df_class_perc_in_sample <- data.frame(
        Event = c(),
        batch = c(),
        mean = c(),
        sd = c()
    )
    for (batch in 1:N_batches) {
        batch_statistics <- simulation_statistics_all[[batch]]
        df_class_perc_in_sample_batch <- data.frame(
            Event = c("Diploid", "WGD", "Missegregation"),
            batch = rep(legend_batches[batch], length = 3),
            mean = c(
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_initial_mean"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_WGD_mean"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_missegregation_mean"), 1]
            ),
            sd = c(
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_initial_sd"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_WGD_sd"), 1],
                batch_statistics$statistics_sample[which(rownames(batch_statistics$statistics_sample) == "perc_missegregation_sd"), 1]
            )
        )
        df_class_perc_in_sample <- rbind(df_class_perc_in_sample, df_class_perc_in_sample_batch)
    }
    #--------------------------Plot the clone counts in total population
    filename <- paste(plotfile_prefix, "_true_clone_count_in_whole_pop.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 500)
    p <- ggplot(df_n_clones_true_in_tot_pop, aes(x = time, y = true_n_clones_whole_pop_mean, group = Parameter, color = Parameter)) +
        geom_line() +
        geom_point() +
        geom_errorbar(
            aes(ymin = true_n_clones_whole_pop_mean - true_n_clones_whole_pop_sd, ymax = true_n_clones_whole_pop_mean + true_n_clones_whole_pop_sd),
            width = .2,
            position = position_dodge(0.05)
        ) +
        coord_cartesian(ylim = c(0, NA))
    # +
    # ylim(0, NA)
    p <- p + labs(
        title = "",
        x = "Age",
        y = "True clone count"
    ) +
        theme_classic()
    print(p)
    dev.off()
    #------------------------------------Plot the clone counts in sample
    filename <- paste(plotfile_prefix, "_clone_count_in_sample.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 500)
    p <- ggplot(df_n_clones_in_sample, aes(x = batch, y = mean, fill = Statistics)) +
        geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
        ) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
            width = .2,
            position = position_dodge(.9)
        ) +
        coord_cartesian(ylim = c(0, NA))
    # +
    # ylim(0, NA)
    p <- p + labs(title = "", x = "Parameter", y = "Count") +
        theme_classic() +
        scale_fill_manual(values = c("#999999", "#E69F00"))
    print(p)
    dev.off()
    #-----------------Plot the percentages of each event class in sample
    filename <- paste(plotfile_prefix, "_event_percentages_in_sample.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 500)
    p <- ggplot(df_class_perc_in_sample, aes(x = batch, y = mean, fill = Event)) +
        geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
        ) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
            width = .2,
            position = position_dodge(.9)
        ) +
        coord_cartesian(ylim = c(0, 100))
    p <- p + labs(title = "", x = "Subgroup", y = "Percentage") +
        theme_classic()
    print(p)
    dev.off()
}
